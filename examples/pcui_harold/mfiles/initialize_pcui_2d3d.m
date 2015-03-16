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

working_folder = '/home/barthur/pcui-3d/examples/pcui_harold';

% -------------------------------------------------------------------------
% Load data from 2d run
% -------------------------------------------------------------------------
load ./2ch_strat_02_2d3d_restart.mat
Nz = 96 + 4*3;

% These are the output files with the filenames stripped out of extensions 
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
fname_rho = 'input_S';
fname_uvw = 'input_UVW';
fname_UVW_to_PCUI = 'uvw_init_from_matlab';
fname_rho_init_to_PCUI = 'rho_init_from_matlab';
fname_rho_full_to_PCUI = 'rho_full_from_matlab';
fname_phi_to_PCUI = 'phi_init_from_matlab';
fname_phi2_to_PCUI = 'phi2_init_from_matlab';
fname_phi3_to_PCUI = 'phi3_init_from_matlab';
fname_phi4_to_PCUI = 'phi4_init_from_matlab';
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

%Reshape to 3d grid and write
%velocity
u_pcui = repmat(u,[1 1 Nz]);
v_pcui = repmat(v,[1 1 Nz]);
w_pcui = repmat(w,[1 1 Nz]);
uvw_pcui = zeros([size(u_pcui),3]);
uvw_pcui(:,:,:,1) = u_pcui;
uvw_pcui(:,:,:,2) = v_pcui;
uvw_pcui(:,:,:,3) = w_pcui;
write_binary_file_pcui(working_folder, fname_UVW_to_PCUI, params, uvw_pcui); 
figure; pcolor(x,y,squeeze(u)); shading flat; colorbar;
figure; pcolor(x,y,squeeze(v)); shading flat; colorbar;
figure; pcolor(x,y,squeeze(w)); shading flat; colorbar;

%density
rho_full_pcui = repmat(rho,[1 1 Nz]); 
write_binary_file_pcui(working_folder, fname_rho_full_to_PCUI, params, rho_full_pcui);
figure; pcolor(x,y,squeeze(rho)); shading flat; colorbar;

%scalars
phi2_pcui = repmat(phi2,[1 1 Nz]); 
write_binary_file_pcui(working_folder, fname_phi2_to_PCUI, params, phi2_pcui);
figure; pcolor(x,y,squeeze(phi2)); shading flat; colorbar;

phi3_pcui = repmat(phi3,[1 1 Nz]); 
write_binary_file_pcui(working_folder, fname_phi3_to_PCUI, params, phi3_pcui);
figure; pcolor(x,y,squeeze(phi3)); shading flat; colorbar;

phi4_pcui = repmat(phi4,[1 1 Nz]); 
write_binary_file_pcui(working_folder, fname_phi4_to_PCUI, params, phi4_pcui);
figure; pcolor(x,y,squeeze(phi4)); shading flat; colorbar;


