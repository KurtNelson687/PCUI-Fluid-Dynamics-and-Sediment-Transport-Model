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

% These are the output files with the filenames stripped out of extensions 
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
fname_rho = 'input_S';
fname_uvw = 'input_UVW';
fname_UVW_to_PCUI = 'uvw_init_from_matlab';
fname_rho_init_to_PCUI = 'rho_init_from_matlab';
fname_rho_full_to_PCUI = 'rho_full_from_matlab';
fname_phi_to_PCUI = 'phi_init_from_matlab';
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

% Load global grid from grid generator file and rearrange y and z
% L = 1; H = 1; W = 1;
% N = 512; 
% dx = L/N; dy = H/N; dz = W/16;
% x = (-1.5*dx:dx:L+1.5*dx)'; x = repmat(x,[1 N+4]); x_global = repmat(x,[1 1 20]);
% y = -H-1.5*dy:dy:1.5*dy;    y = repmat(y,[N+4 1]); y_global = repmat(y,[1 1 20]);
% z = -1.5*dz:dz:W+1.5*dz;    z = reshape(z,[1 1 20]); z_global = repmat(z,[N+4 N+4 1]);
% load '/home/barthur/zang/grids/grid_energy_test.mat'
load '/home/barthur/zang/grids/grid_1152x128x64_r102_w1_zstretch_s218_f128.mat';
% load '/home/barthur/zang/grids/grid_1024x128x128_r102_w125_zstretch.mat'
% load '/home/barthur/zang/grids/grid_2D_test.mat'
% load '/home/barthur/zang/grids/grid_local_2D_test.mat';
x_global = x; x_global = permute(x_global,[1 3 2]);  %DONT FORGET THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
y_global = z; y_global = permute(y_global,[1 3 2]);
z_global = y; z_global = permute(z_global,[1 3 2]);

%Prepare for writing (account for multiple processors with 2 cell halos)
ni_wh = params.ni+4*params.px; %includes halo for each proc
nj_wh = params.nj+4*params.py;
nk_wh = params.nk+4*params.pz;
x_pcui = zeros(ni_wh,nj_wh,nk_wh);
y_pcui = x_pcui; z_pcui = x_pcui;

nni = params.ni/params.px+4;
nnj = params.nj/params.py+4;
nnk = params.nk/params.pz+4;
for idk = 0:params.pz-1
    for idj = 0:params.py-1
        for idi = 0:params.px-1 
            indx_i_pcui_1     = idi*nni+1;
            indx_i_pcui_end   = (idi+1)*nni;
            indx_j_pcui_1     = idj*nnj+1;
            indx_j_pcui_end   = (idj+1)*nnj;
            indx_k_pcui_1     = idk*nnk+1;
            indx_k_pcui_end   = (idk+1)*nnk;
            
            indx_i_global_1   = indx_i_pcui_1 - idi*4;
            indx_i_global_end = indx_i_global_1 + nni - 1;
            indx_j_global_1   = indx_j_pcui_1 - idj*4;
            indx_j_global_end = indx_j_global_1 + nnj - 1;
            indx_k_global_1   = indx_k_pcui_1 - idk*4;
            indx_k_global_end = indx_k_global_1 + nnk - 1;
 
            x_pcui(indx_i_pcui_1:indx_i_pcui_end, ...
                   indx_j_pcui_1:indx_j_pcui_end, ...
                   indx_k_pcui_1:indx_k_pcui_end) = ...
            x_global(indx_i_global_1:indx_i_global_end, ...
                   indx_j_global_1:indx_j_global_end, ...
                   indx_k_global_1:indx_k_global_end);
            y_pcui(indx_i_pcui_1:indx_i_pcui_end, ...
                   indx_j_pcui_1:indx_j_pcui_end, ...
                   indx_k_pcui_1:indx_k_pcui_end) = ...
            y_global(indx_i_global_1:indx_i_global_end, ...
                   indx_j_global_1:indx_j_global_end, ...
                   indx_k_global_1:indx_k_global_end);
            z_pcui(indx_i_pcui_1:indx_i_pcui_end, ...
                   indx_j_pcui_1:indx_j_pcui_end, ...
                   indx_k_pcui_1:indx_k_pcui_end) = ...
            z_global(indx_i_global_1:indx_i_global_end, ...
                   indx_j_global_1:indx_j_global_end, ...
                   indx_k_global_1:indx_k_global_end);    
        end
    end
end

xyz_pcui = zeros(ni_wh,nj_wh,nk_wh,3);
xyz_pcui(:,:,:,1) = x_pcui;
xyz_pcui(:,:,:,2) = y_pcui;
xyz_pcui(:,:,:,3) = z_pcui;

% Write PCUI binary files depending on the number of processors
write_binary_file_pcui(working_folder, fname_grid_to_PCUI, params, xyz_pcui);

% -------------------------------------------------------------------------
% Initialize PCUI with a solitary wave
% -------------------------------------------------------------------------
% Prepare density field
% h1 = -0.3;
% a = 0.05;
% Lw = 0.9;
% delta = 0.02;
% alpha = 0.99;
% rho_init_pcui = ones(size(x_pcui));
% zeta = -a*exp(-(x_pcui/Lw).^2); % + 0.001*rand(size(x_pcui));
% % zeta = -a*sech(x_pcui/Lw).^2 + 0.001*rand(size(x_pcui));
% rho_pert_pcui = -0.5*0.03*tanh(2*(y_pcui - zeta - h1)/delta*atanh(alpha));
% rho_full_pcui = rho_init_pcui+rho_pert_pcui;
% u_pcui = zeros(size(rho_init_pcui));
% v_pcui = u_pcui; w_pcui = u_pcui;
% uvw_pcui(:,:,:,1) = u_pcui;
% uvw_pcui(:,:,:,2) = v_pcui;
% uvw_pcui(:,:,:,3) = w_pcui;

% Write PCUI binary files depending on the number of processors
% write_binary_file_pcui(working_folder, fname_rho_full_to_PCUI, params, rho_full_pcui);
% write_binary_file_pcui(working_folder, fname_rho_init_to_PCUI, params, rho_init_pcui);
% write_binary_file_pcui(working_folder, fname_UVW_to_PCUI, params, uvw_pcui); 

% -------------------------------------------------------------------------
% Initialize PCUI with a passive scalar
% -------------------------------------------------------------------------
% phi_pcui = zeros(size(x_pcui));
% phi_pcui(x_pcui>3.52 & x_pcui<4) = 1;
% phi_pcui = x_pcui;
% write_binary_file_pcui(working_folder, fname_phi_to_PCUI, params, phi_pcui);

% -------------------------------------------------------------------------
% Verify initialized solitary wave
% -------------------------------------------------------------------------
% Plot density field
% x_plot = squeeze(x_global(3:end-2,3:end-2,1));
% y_plot = squeeze(y_global(3:end-2,3:end-2,1));
% 
% fig1 = figure(1);
% clf
% set(fig1,'Renderer','zbuffer');
% set(fig1,'Color','white');
% rho_init_plot = ones(size(x_plot));
% zeta_plot = -a*exp(-(x_plot/Lw).^2) + 0.001*rand(size(x_plot));
% % zeta_plot = -a*sech(x_plot/Lw).^2 + 0.001*rand(size(x_plot));
% rho_pert_plot = -0.5*0.03*tanh(2*(y_plot - zeta_plot - h1)/delta*atanh(alpha));
% rho_full_plot = rho_init_plot+rho_pert_plot;
% pcolor(x_plot,y_plot,rho_full_plot);
% axis image;
% shading flat;
% colorbar;

% -------------------------------------------------------------------------
% Verify initialized passive scalar 
% -------------------------------------------------------------------------
% phi_plot = zeros(size(x_plot));
% phi_plot(x_plot>3.52 & x_plot<4) = 1;
% phi_plot = x_plot;
% 
% fig2 = figure(2);
% clf
% set(fig2,'Renderer','zbuffer');
% set(fig2,'Color','white');
% pcolor(x_plot,y_plot,phi_plot);
% axis image;
% shading flat;
% colorbar;

% -------------------------------------------------------------------------
% Verify initialized grid
% -------------------------------------------------------------------------
%Plot grid
fig3 = figure(3);
clf
set(fig3,'Renderer','zbuffer');
set(fig3,'Color','white');
plot(squeeze(x_global(:,:,1)),squeeze(y_global(:,:,1)),'k.');
xlabel('x [m]');
ylabel('y [m]');
axis equal;

fig4 = figure(4);
clf
set(fig4,'Renderer','zbuffer');
set(fig4,'Color','white');
plot(squeeze(x_global(:,1,:)),squeeze(z_global(:,1,:)),'k.');
xlabel('x [m]');
ylabel('z [m]');
axis equal;

fig5 = figure(5);
clf
set(fig5,'Renderer','zbuffer');
set(fig5,'Color','white');
plot(squeeze(z_global(1,:,:)),squeeze(y_global(1,:,:)),'k.');
xlabel('z [m]');
ylabel('y [m]');
axis equal;