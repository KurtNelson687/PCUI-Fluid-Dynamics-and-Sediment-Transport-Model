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

working_folder = '../';

% These are the output files with the filenames stripped out of extensions 
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
fname_rho = 'input_S';
fname_uvw = 'input_UVW';
fname_UVW_to_PCUI = 'uvw_init_from_matlab';
fname_rho_init_to_PCUI = 'rho_init_from_matlab';
fname_rho_full_to_PCUI = 'rho_full_from_matlab';
fname_u_west_to_pcui = 'u_west_init_from_matlab';

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
uvw_pcui = zeros([size(x),3]);
rho_init_pcui = ones(size(x));
rho_full_pcui = rho_init_pcui;
% u_prof = 0.1*ones(size(y(1,:,1)')); %uniform flow
% u_prof = 0.05/0.3*y(1,:,1)'; %pure shear flow
load fastfit %load fast data
u_prof = [ftanh(y(1,y(1,:,1)<jhoverlap,1)); flog(y(1,y(1,:,1)>=jhoverlap,1))];
u_west_pcui = repmat(u_prof,[1 params.ni+4*params.px params.nk+4*params.pz]);
u_west_pcui = permute(u_west_pcui,[2,1,3]);
u_pcui = repmat(u_prof,[1 params.ni+4*params.px params.nk+4*params.pz]);
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
write_binary_file_pcui(working_folder, fname_u_west_to_pcui, params, u_west_pcui);

% fid = fopen('../u_west_init_from_matlab.700', 'w');
% precision = 'float64';
% headerflag = 'float32'; % This needs to be changed to float64 to read in C binary data
% footerflag = headerflag;
% a = fwrite(fid,1,headerflag);%header
% data = fwrite(fid,u_west_pcui,precision);%data
% b = fwrite(fid,1, footerflag);%footer

% % -------------------------------------------------------------------------
% % Verify initialized density
% % -------------------------------------------------------------------------
% fig1 = figure(1);
% clf
% set(fig1,'Renderer','zbuffer');
% set(fig1,'Color','white');
% rho_init_plot = ones(size(x_plot));
% rho_full_plot = rho_init_plot;
% pcolor(x_plot,y_plot,rho_full_plot);
% colorbar;
% colormap gray;
% axis equal;
% 
% -------------------------------------------------------------------------
% Verify initialized velocity profile (same as u_west)
% -------------------------------------------------------------------------
fig2 = figure(2);
clf
set(fig2,'Renderer','zbuffer');
set(fig2,'Color','white');
% u_prof_plot = 0.1*ones(size(y_plot(1,:,1)'));
% u_prof_plot = 0.05/0.3*y_plot(1,:,1)';
u_prof_plot = [ftanh(y_plot(1,y_plot(1,:,1)<jhoverlap,1)); flog(y_plot(1,y_plot(1,:,1)>=jhoverlap,1))];
plot(u_prof_plot,y_plot(1,:),'k-o');
axis([0 max(u_prof_plot) 0 D]);

% -------------------------------------------------------------------------
% Verify initialized initial velocity
% -------------------------------------------------------------------------
fig3 = figure(3);
clf
set(fig3,'Renderer','zbuffer');
set(fig3,'Color','white');
% u_plot = 0.1*ones(size(y_plot(1,:,1)'));
% u_plot = 0.05*3*y_plot(1,:,1)';
u_plot = [ftanh(y_plot(1,y_plot(1,:,1)<jhoverlap,1)); flog(y_plot(1,y_plot(1,:,1)>=jhoverlap,1))];
u_plot = repmat(u_plot,[1 params.ni]);
u_plot = permute(u_plot,[2,1]);
pcolor(x_plot,y_plot,u_plot);
colorbar;
shading flat;
colormap gray;
axis equal;
