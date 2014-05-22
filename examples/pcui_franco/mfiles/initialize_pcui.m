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

working_folder = '/Users/fjzarama/Desktop/Research/DownstreamDevelopment/pcui-3d/examples/pcui_franco';

% These are the output files with the filenames stripped out of extensions 
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
fname_rho = 'input_S';
fname_uvw = 'input_UVW';
fname_UVW_to_PCUI = 'uvw_init_from_matlab';
fname_rho_init_to_PCUI = 'rho_init_from_matlab';
fname_rho_full_to_PCUI = 'rho_full_from_matlab';

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
% (includes ghost grid points of each submatrix)
[x,y,z] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,1);
% read the files containing the grid definition and assemble into array
% (does not include ghost grid points of each submatrix)
[x_plot,y_plot,z_plot] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,0);
x_plot = squeeze(x_plot(:,:,floor(length(z_plot(1,1,:)/2))));
y_plot = squeeze(y_plot(:,:,floor(length(z_plot(1,1,:)/2))));

% -------------------------------------------------------------------------
% Initialize PCUI with constant density
% -------------------------------------------------------------------------
% Prepare density field
D = params.by;
rho_init_pcui = ones(size(x));
rho_full_pcui = rho_init_pcui;
% u_pcui = 0.05/0.3*y; %pure shear flow
load fastfit %load fast data
u_pcui = [ftanh(y(1,y(1,:,1)<jhoverlap,1)); flog(y(1,y(1,:,1)>=jhoverlap,1))];
u_pcui = repmat(u_pcui,[1 144 20]);
u_pcui = permute(u_pcui,[2,1,3]);
v_pcui = zeros(size(rho_init_pcui)); 
w_pcui = zeros(size(rho_init_pcui));
uvw_pcui(:,:,:,1) = u_pcui;
uvw_pcui(:,:,:,2) = v_pcui;
uvw_pcui(:,:,:,3) = w_pcui;

% Write PCUI binary files depending on the number of processors
write_binary_file_pcui(working_folder, fname_rho_full_to_PCUI, params, rho_full_pcui);
write_binary_file_pcui(working_folder, fname_rho_init_to_PCUI, params, rho_init_pcui);
write_binary_file_pcui(working_folder, fname_UVW_to_PCUI, params, uvw_pcui); 

% -------------------------------------------------------------------------
% Verify initialized density
% -------------------------------------------------------------------------
% Plot density field
fig1 = figure(1);
clf
set(fig1,'Renderer','zbuffer');
set(fig1,'Color','white');
rho_init_plot = ones(size(x_plot));
rho_full_plot = rho_init_plot;
pcolor(x_plot,y_plot,rho_full_plot);
colorbar;
axis equal;
