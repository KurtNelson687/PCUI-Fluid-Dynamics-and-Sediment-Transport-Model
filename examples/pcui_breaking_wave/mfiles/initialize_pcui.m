%
% Filename    : initialize_pcui.m
% Author      : Goncalo Gil
% Description : Creates an input file to be read by PCUI. It can write
%               three- or four-dimensional arrays for the density and 
%               velocity fields, respectively.                           
%
% Author      : Goncalo Gil, Stanford University
% email       : gilg@stanford.edu
%

clear all; close all; clc;

working_folder = '/home/barthur/pcui-3d/examples/pcui_breaking_wave';

% These are the output files with the filenames stripped out of extensions 
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
fname_rho = 'input_S';
fname_uvw = 'input_UVW';
fname_UVW_to_PCUI = 'uvw_init_from_matlab';
fname_rho_init_to_PCUI = 'rho_init_from_matlab';
fname_rho_full_to_PCUI = 'rho_full_from_matlab';
fname_grid_to_PCUI = 'xyz_init_from_matlab';

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

% -------------------------------------------------------------------------
% Initialize PCUI grid
% -------------------------------------------------------------------------
load /home/barthur/zang/grids/pcui_test.mat
x_pcui_global = x; x_pcui_global = permute(x_pcui_global,[1 3 2]);
y_pcui_global = z; y_pcui_global = permute(y_pcui_global,[1 3 2]);
z_pcui_global = y; z_pcui_global = permute(z_pcui_global,[1 3 2]);

x_pcui = x_pcui_global;
y_pcui = y_pcui_global;
z_pcui = z_pcui_global;

xyz_pcui = zeros(params.ni+4,params.nj+4,params.nk+4,3);
xyz_pcui(:,:,:,1) = x_pcui;
xyz_pcui(:,:,:,2) = y_pcui;
xyz_pcui(:,:,:,3) = z_pcui;

% Write PCUI binary files depending on the number of processors
write_binary_file_pcui(working_folder, fname_grid_to_PCUI, params, xyz_pcui);

% -------------------------------------------------------------------------
% Initialize PCUI with a solitary wave
% -------------------------------------------------------------------------
% Prepare density field
h1 = -0.3;
a = 0.1;
Lw = 0.7;
delta = 0.2;
alpha = 0.99;
rho_init_pcui = ones(size(x_pcui));
zeta = -a*exp(-(x_pcui/Lw).^2); % + 0.001*randn(size(x));
rho_pert_pcui = -0.5*0.03*tanh(2*(y_pcui - zeta - h1)/delta*atanh(alpha));
rho_full_pcui = rho_init_pcui+rho_pert_pcui;
u_pcui = zeros(size(rho_init_pcui));
v_pcui = u_pcui; w_pcui = u_pcui;
uvw_pcui(:,:,:,1) = u_pcui;
uvw_pcui(:,:,:,2) = v_pcui;
uvw_pcui(:,:,:,3) = w_pcui;

% Write PCUI binary files depending on the number of processors
write_binary_file_pcui(working_folder, fname_rho_full_to_PCUI, params, rho_full_pcui);
write_binary_file_pcui(working_folder, fname_rho_init_to_PCUI, params, rho_init_pcui);
write_binary_file_pcui(working_folder, fname_UVW_to_PCUI, params, uvw_pcui); 

% -------------------------------------------------------------------------
% Verify initialized solitary wave
% -------------------------------------------------------------------------
% Plot density field
x_plot = squeeze(x_pcui_global(3:end-2,3:end-2,1));
y_plot = squeeze(y_pcui_global(3:end-2,3:end-2,1));

fig1 = figure(1);
clf
set(fig1,'Renderer','zbuffer');
set(fig1,'Color','white');
rho_init_plot = ones(size(x_plot));
zeta_plot = -a*exp(-(x_plot/Lw).^2); % + 0.001*randn(size(x_plot));
rho_pert_plot = -0.5*0.03*tanh(2*(y_plot - zeta_plot - h1)/delta*atanh(alpha));
rho_full_plot = rho_init_plot+rho_pert_plot;
pcolor(x_plot,y_plot,rho_full_plot);
axis image;
colorbar;

% -------------------------------------------------------------------------
% Verify initialized grid
% -------------------------------------------------------------------------
%Plot grid
fig2 = figure(2);
clf
set(fig2,'Renderer','zbuffer');
set(fig2,'Color','white');
plot(squeeze(x_pcui(:,:,1)),squeeze(y_pcui(:,:,1)),'k.');
xlabel('x [m]');
ylabel('y [m]');
axis equal;

fig3 = figure(3);
clf
set(fig3,'Renderer','zbuffer');
set(fig3,'Color','white');
plot(squeeze(x_pcui(:,1,:)),squeeze(z_pcui(:,1,:)),'k.');
xlabel('x [m]');
ylabel('z [m]');
axis equal;

fig4 = figure(4);
clf
set(fig4,'Renderer','zbuffer');
set(fig4,'Color','white');
plot(squeeze(z_pcui(1,:,:)),squeeze(y_pcui(1,:,:)),'k.');
xlabel('z [m]');
ylabel('y [m]');
axis equal;