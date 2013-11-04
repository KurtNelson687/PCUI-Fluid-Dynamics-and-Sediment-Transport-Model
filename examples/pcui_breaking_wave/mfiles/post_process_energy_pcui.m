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

smooth = 0; % Interpolate shading in plot
isprint = 0; % Print plot to eps file
plot_name = 'breaking_wave';
working_folder = '../';
save_folder =    '../figs';
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
params.eddy_viscosity = variable_value_pcui('ak',ftext);
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
[Ep,Eb,Ea,dEpdt,dEbdt,dEadt,phi_d] = read_binary_energy_pcui(working_folder, fname_energy, params.dt);

% -------------------------------------------------------------------------
% Calculate dissipation
% -------------------------------------------------------------------------
iskip = 1;
istart =  1;
iend = floor(params.nsteps/params.nsave);

eps = zeros(1,iend); Ek = eps; phi_z = eps;
eps_top = eps; eps_int = eps; eps_low = eps; eps_bot = eps; eps_tot = eps;
for istep = istart:iskip:iend
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
    
    %Calculate phi_z
    phi_z(istep) = 9.81*sum(sum(sum(rho.*v(2:end-1,2:end-1,2:end-1).*metrics.J)));

    %Calculate dissipation
%     eps(istep) = calculate_binary_dissipation(u,v,w,10^-6,metrics);
    [eps_top(istep),eps_int(istep),eps_low(istep),eps_bot(istep),eps_tot(istep)] = ...
               calculate_binary_dissipation_breakdown(u,v,w,rho,10^-6,metrics);
    
    %Calculate kinetic energy
    Ek(istep) = 0.5*sum(sum(sum((u(2:end-1,2:end-1,2:end-1).^2 ...
                                +v(2:end-1,2:end-1,2:end-1).^2 ...
                                +w(2:end-1,2:end-1,2:end-1).^2).*metrics.J)));
end

Et = Ek + Ep(1:params.nsave:params.nsteps)';
dEtdt = (Et(3:end)-Et(1:end-2))/2/params.dt/params.nsave;
dEkdt = (Ek(3:end)-Ek(1:end-2))/2/params.dt/params.nsave;

%%
% -------------------------------------------------------------------------
% Plot energy quantities
% -------------------------------------------------------------------------
t = params.dt:params.dt:params.dt*params.nsteps;
tdt = t(2:end-1);
tn = params.dt*istart:params.dt*params.nsave:params.dt*params.nsteps;

figure(1);
plot(tn,Et,'k',t,Ep,'k--',t,Eb,'b',t,Ea,'r',tn,Ek,'k:');

figure(2);
plot(tn(2:end-1),dEtdt,'k',tdt,dEpdt,'k--',tdt,dEbdt,'b',tdt,dEadt,'r',tn(2:end-1),dEkdt,'k:',tn,phi_z,'m')

figure(3);
% plot(t,phi_d,'b',tn,eps,'k',tdt,dEbdt*10^-1,'b');
plot(tn,eps_top,'m',tn,eps_int,'b',tn,eps_low,'g',tn,eps_bot,'r',tn,eps_tot,'k')
% bar(tn,[eps_bot' eps_int' eps_top' eps_low'],'stack')
