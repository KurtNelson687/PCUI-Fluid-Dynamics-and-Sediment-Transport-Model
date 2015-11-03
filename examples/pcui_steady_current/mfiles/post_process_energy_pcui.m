%
% Filename    : post_process_pcui.m
% Author      : Goncalo Gil
% Description : Reads in PCUI data and plots a slice of the velocity field
%               and density field for the internal seiche test run.
%                                       
% Author      : Goncalo Gil, Stanford University
% email       : gilg@stanford.edu
%

clear all; close all; clc;

working_folder = './..';
save_file = 'energy_256x256_k0';
TEND = 1200;

% These are the output files with the filenames stripped out of extensions 
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
fname_rho = 'output_S';
fname_uvw = 'output_UVW';
fname_energy = 'output_energy';

% -------------------------------------------------------------------------
% Get problem parameters and variables from PCUI
% -------------------------------------------------------------------------
% read the file containing the parameter definition
ftext = fileread(fullfile(working_folder, 'io.f'));
params.dt = variable_value_pcui('dtime',ftext);
params.molecular_viscosity = variable_value_pcui('vis',ftext);
params.molecular_diffusivity = variable_value_pcui('ak',ftext);
params.nsteps = variable_value_pcui('nstep',ftext);
params.nsave = variable_value_pcui('nsave',ftext);

% read the file containing the domain definition
ftext = fileread(fullfile(working_folder, 'cavity.f'));
params.bx = variable_value_pcui('bx',ftext);
params.by = variable_value_pcui('by',ftext);
params.bz = variable_value_pcui('bz',ftext);

% read the file containing the grid size and processor definitions
ftext = fileread(fullfile(working_folder, 'size.inc'));
params.ni = variable_value_pcui('ni',ftext);
params.nj = variable_value_pcui('nj',ftext);
params.nk = variable_value_pcui('nk',ftext);
params.px = variable_value_pcui('px',ftext);
params.py = variable_value_pcui('py',ftext);
params.pz = variable_value_pcui('pz',ftext);

% read the files containing the grid definition and assemble into array
% (includes 2 halo cells for each proc)
[x,y,z] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,1);

%extract physical grid points including 1 halo cell
nni = params.ni/params.px+4;
nnj = params.nj/params.py+4;
nnk = params.nk/params.pz+4;
exi = []; exj = []; exk = [];
if params.px == 1
    exi = 2:nni-1;
else
    for i=0:params.px-1
        if i==0
            exi = [exi,2:nni-2];
        elseif i==params.px-1
            exi = [exi,i*nni+3:(i+1)*nni-1];
        else
            exi = [exi,i*nni+3:(i+1)*nni-2];
        end
    end
end

if params.py == 1
    exj = 2:nnj-1;
else
    for j=0:params.py-1
        if j==0
            exj = [exj,2:nnj-2];
        elseif j==params.py-1
            exj = [exj,j*nnj+3:(j+1)*nnj-1];
        else
            exj = [exj,j*nnj+3:(j+1)*nnj-2];
        end
    end
end

if params.pz == 1
    exk = 2:nnk-1;
else
    for k=0:params.pz-1
        if k==0
            exk = [exk,2:nnk-2];
        elseif k==params.pz-1
            exk = [exk,k*nnk+3:(k+1)*nnk-1];
        else
            exk = [exk,k*nnk+3:(k+1)*nnk-2];
        end
    end
end
    
x = x(exi,exj,exk);
y = y(exi,exj,exk);
z = z(exi,exj,exk);

%calculate metric quantities
metrics = calculate_binary_metrics(x,y,z);

% read the energy file
[Ep,Eb,Ea,dEpdt,dEbdt,dEadt,phi_d] = read_binary_energy_pcui( ...
                                  working_folder, fname_energy, params.dt);

% -------------------------------------------------------------------------
% Calculate other energy quantities
% -------------------------------------------------------------------------
istart = 2; %after initial condition output
%iend = params.nsteps/params.nsave + 1;
iend = TEND/params.nsave + 1;

nu = params.molecular_viscosity; 
kappa = params.molecular_diffusivity;
g = 9.81;

%Ek = zeros(params.nsteps/params.nsave,1); 
Ek = zeros(TEND/params.nsave,1); 
phi_z = Ek; phi_i = Ek;
%eps_top = Ek; eps_int = Ek; eps_low = Ek; eps_bot = Ek; eps_tot = Ek; eps_both = Ek;
for istep = istart:iend
    display(istep);
    
    % Read velocity field (includes 1 halo cell)
    [u,v,w] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
                                      params, 1,1);     
    u = u(exi,exj,exk);
    v = v(exi,exj,exk);
    w = w(exi,exj,exk);
    
    % Read density field (includes 0 halo cells)
    rho = read_binary_file_pcui(working_folder, fname_rho, istep, ...
                                 params, 0,0);   
                             
    %Calculate phi_i
    phi_i(istep-1) = -kappa*g*...
            ( sum(sum(squeeze(rho(:,end,:)).*squeeze(metrics.ET_Y(:,end,:)))) ...
             -sum(sum(squeeze(rho(:,1,  :)).*squeeze(metrics.ET_Y(:,1,  :)))) );
    
    %Calculate phi_z
    phi_z(istep-1) = g*sum(sum(sum(rho.*v(2:end-1,2:end-1,2:end-1).*metrics.J)));

%     %Calculate dissipation
%     [eps_top(istep-1),eps_int(istep-1),eps_low(istep-1),eps_bot(istep-1),eps_tot(istep-1),eps_both(istep-1)] = ...
%                calculate_binary_dissipation_breakdown(u,v,w,rho,delta_rho,nu,delta_s,y(2:end-1,2:end-1,2:end-1),metrics);
%    
    %Calculate kinetic energy
    Ek(istep-1) = 0.5*sum(sum(sum((u(2:end-1,2:end-1,2:end-1).^2 ...
                                +v(2:end-1,2:end-1,2:end-1).^2 ...
                                +w(2:end-1,2:end-1,2:end-1).^2).*metrics.J)));
end

Et = Ek + Ep(params.nsave:params.nsave:TEND);
dEtdt = (Et(3:end)-Et(1:end-2))/2/params.dt/params.nsave;
dEkdt = (Ek(3:end)-Ek(1:end-2))/2/params.dt/params.nsave;
phi_e = phi_d(params.nsave:params.nsave:TEND) - phi_i;

% -------------------------------------------------------------------------
% Save energy quantities
% -------------------------------------------------------------------------
t = params.dt*(1:TEND);
tdt = t(2:end-1);
tn = params.dt*(params.nsave:params.nsave:TEND);
nsave = params.nsave;
nsteps = TEND;
save(save_file,'t','tdt','tn','nsave','nsteps','Ek','Ep','Eb','Ea', ...
               'dEpdt','dEbdt','dEadt','dEkdt','dEtdt',...
               'phi_d','phi_z','phi_i','phi_e');
display('Complete');

