% Displays PCUI particle tracking data 
clear all; clc; close all;

% PLOTTING OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
working_folder = '../';
% These are the output files with the filenames stripped out of extensions 
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
fname_rho = 'output_S';
fname_phi = 'output_phi';
fname_uvw = 'output_UVW';
filename_xpart = 'output_xPart.dat';
filename_upart = 'output_uPart.dat';
fname_vst = 'output_vst_o';
fname_akst = 'output_akst_o';
fname_diss_sgs = 'output_diss_sgs';

% OTHER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the file containing the parameter definition
ftext = fileread(fullfile(working_folder, 'io.f'));
params.dt = variable_value_pcui('dtime',ftext);
params.molecular_viscosity = variable_value_pcui('vis',ftext);
params.eddy_viscosity = variable_value_pcui('ak',ftext);
params.nsteps = variable_value_pcui('nstep',ftext);
params.nsave = variable_value_pcui('nsave',ftext);

% read the file containing the domain definition
ftext = fileread(fullfile(working_folder, 'cavity.f'));
params.bx = variable_value_pcui('bx',ftext); x_length = params.bx;
params.by = variable_value_pcui('by',ftext); z_length = params.by;
params.bz = variable_value_pcui('bz',ftext); y_length = params.bz;

% read the file containing the grid size and processor definitions
ftext = fileread(fullfile(working_folder, 'size.inc'));
params.ni = variable_value_pcui('ni',ftext); Nx = params.ni;
params.nj = variable_value_pcui('nj',ftext); Nz = params.nj;
params.nk = variable_value_pcui('nk',ftext); Ny = params.nk;
params.px = variable_value_pcui('px',ftext);
params.py = variable_value_pcui('py',ftext);
params.pz = variable_value_pcui('pz',ftext);

%Load grid
[x,y,z] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,0);
x = squeeze(x(:,:,1));
y = squeeze(y(:,:,1));

%Find correct istep value
n = [0, params.nsave:params.nsave:params.nsteps, params.nsteps+1];

TEND = 400;
xpartall = nan(params.ni*params.nj,3,1);
for timestep = 0:params.nsave:TEND
    display(timestep);
    istep = find(n==timestep,1);
    xpartall(:,:,istep) = read_binary_particles_pcui(working_folder, filename_xpart, istep, params);
end

%%
clc; close all;
np = params.ni;
figure;
hold on;
for timestep=0:params.nsave:TEND
    istep = find(n==timestep,1);
    display(istep);
    rho = read_binary_file_pcui(working_folder, fname_rho, istep, ...
                                 params, 0,0);  
    cla;
    contour(x,y,squeeze(rho(:,:,1)),[1 1],'r');
    plot(squeeze(xpartall(:,1,1:istep)),squeeze(xpartall(:,2,1:istep)),'k.');
    axis equal;
    axis([0 params.bx 0 params.by]);
    drawnow;
    pause;
end
hold off;
    