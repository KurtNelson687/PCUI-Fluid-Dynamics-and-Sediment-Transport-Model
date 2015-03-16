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
% load './grid_case3_2D.mat'
% load '/home/barthur/zang/grids/grid_energy_test.mat'
% load './2ch_flat_part.mat'
load '/home/barthur/zang/grids/grid_1536x128x16_r102_w1_zstretch_s218_f128_djl.mat';
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

%%
% -------------------------------------------------------------------------
% Initialize PCUI with a solitary wave
% -------------------------------------------------------------------------
% Prepare background density field
load /home/barthur/software/djl-public/djl_new.mat
rho_min = 985;
rho_0 = 1000;
rho_max = 1015;
gamma = 150;
y0 = 0.45;
H = 0.75;
rho_full_pcui = rho_min + 0.5*(rho_max-rho_min)*(1+tanh(gamma*(y0-(y_pcui+H))));

%Interpolate high-res DJL solution to lower-res grid and place in rho_full_pcui
N_flat = 128;
procs_flat = 4;

xdjltmp = xdjl(:);
ydjltmp = ydjl(:) - H;
rhodjltmp = rhodjl(:);
xpcuitmp = x_pcui(1:N_flat+procs_flat*4,:,:); xpcuitmp = xpcuitmp(:);
ypcuitmp = y_pcui(1:N_flat+procs_flat*4,:,:); ypcuitmp = ypcuitmp(:);
rhodjl_low = griddata(xdjltmp,ydjltmp,rhodjltmp,xpcuitmp,ypcuitmp,'cubic');
rho_full_pcui(1:N_flat+procs_flat*4,:,:) = reshape(rhodjl_low,[N_flat+procs_flat*4,nj_wh,nk_wh]);

%adjust BCs
rho_full_pcui(1:2,:,:) = repmat(rho_full_pcui(3,:,:),[2 1 1]); %left
rho_full_pcui(N_flat+procs_flat*4-1:N_flat+procs_flat*4,:,:) = repmat(rho_full_pcui(N_flat+procs_flat*4-2,:,:),[2 1 1]); %right
rho_full_pcui(1:N_flat+procs_flat*4,1:2,:) = repmat(rho_full_pcui(1:N_flat+procs_flat*4,3,:),[1 2 1]); %bottom
rho_full_pcui(1:N_flat+procs_flat*4,end-1:end,:) = repmat(rho_full_pcui(1:N_flat+procs_flat*4,end-2,:),[1 2 1]); %top
rho_full_pcui(1:N_flat+procs_flat*4,:,1:2) = repmat(rho_full_pcui(1:N_flat+procs_flat*4,:,3),[1 1 2]); %front
rho_full_pcui(1:N_flat+procs_flat*4,:,end-1:end) = repmat(rho_full_pcui(1:N_flat+procs_flat*4,:,end-2),[1 1 2]); %back

%normalize
rho_full_pcui = rho_full_pcui/rho_0;

%plot
figure;
pcolor(squeeze(x_pcui(:,:,1)),squeeze(y_pcui(:,:,1)),squeeze(rho_full_pcui(:,:,1)));
shading flat;
colorbar;
% for k=1:nk_wh
%     display(k);
%     cla;
%     pcolor(squeeze(x_pcui(:,:,k)),squeeze(y_pcui(:,:,k)),squeeze(rho_full_pcui(:,:,k)));
%     shading flat;
%     colorbar; 
%     drawnow;
%     pause;
% end

% Prepare velocity field
u_pcui = zeros(size(rho_full_pcui));
v_pcui = u_pcui; w_pcui = u_pcui;

%Interpolate high-res DJL solution to lower-res grid and place in u_pcui/vpcui
udjltmp = udjl(:);
vdjltmp =vdjl(:);
udjl_low = griddata(xdjltmp,ydjltmp,udjltmp,xpcuitmp,ypcuitmp,'cubic');
vdjl_low = griddata(xdjltmp,ydjltmp,vdjltmp,xpcuitmp,ypcuitmp,'cubic');
u_pcui(1:N_flat+procs_flat*4,:,:) = reshape(udjl_low,[N_flat+procs_flat*4,nj_wh,nk_wh]);
v_pcui(1:N_flat+procs_flat*4,:,:) = reshape(vdjl_low,[N_flat+procs_flat*4,nj_wh,nk_wh]);

%adjust BCs - u
u_pcui(1:2,:,:) = repmat(u_pcui(3,:,:),[2 1 1]); %left
u_pcui(N_flat+procs_flat*4-1:N_flat+procs_flat*4,:,:) = repmat(u_pcui(N_flat+procs_flat*4-2,:,:),[2 1 1]); %right
u_pcui(1:N_flat+procs_flat*4,1:2,:) = repmat(u_pcui(1:N_flat+procs_flat*4,3,:),[1 2 1]); %bottom
u_pcui(1:N_flat+procs_flat*4,end-1:end,:) = repmat(u_pcui(1:N_flat+procs_flat*4,end-2,:),[1 2 1]); %top
u_pcui(1:N_flat+procs_flat*4,:,1:2) = repmat(u_pcui(1:N_flat+procs_flat*4,:,3),[1 1 2]); %front
u_pcui(1:N_flat+procs_flat*4,:,end-1:end) = repmat(u_pcui(1:N_flat+procs_flat*4,:,end-2),[1 1 2]); %back

%adjust BCs -v
v_pcui(1:2,:,:) = repmat(v_pcui(3,:,:),[2 1 1]); %left
v_pcui(N_flat+procs_flat*4-1:N_flat+procs_flat*4,:,:) = repmat(v_pcui(N_flat+procs_flat*4-2,:,:),[2 1 1]); %right
v_pcui(1:N_flat+procs_flat*4,1:2,:) = repmat(v_pcui(1:N_flat+procs_flat*4,3,:),[1 2 1]); %bottom
v_pcui(1:N_flat+procs_flat*4,end-1:end,:) = repmat(v_pcui(1:N_flat+procs_flat*4,end-2,:),[1 2 1]); %top
v_pcui(1:N_flat+procs_flat*4,:,1:2) = repmat(v_pcui(1:N_flat+procs_flat*4,:,3),[1 1 2]); %front
v_pcui(1:N_flat+procs_flat*4,:,end-1:end) = repmat(v_pcui(1:N_flat+procs_flat*4,:,end-2),[1 1 2]); %back

%plot
figure;
pcolor(squeeze(x_pcui(:,:,1)),squeeze(y_pcui(:,:,1)),squeeze(u_pcui(:,:,1)));
shading flat;
colorbar;
figure;
pcolor(squeeze(x_pcui(:,:,1)),squeeze(y_pcui(:,:,1)),squeeze(v_pcui(:,:,1)));
shading flat;
colorbar;
% for k=1:nk_wh
%     display(k);
%     cla;
%     pcolor(squeeze(x_pcui(:,:,k)),squeeze(y_pcui(:,:,k)),squeeze(w_pcui(:,:,k)));
%     shading flat;
%     colorbar; 
%     drawnow;
%     pause;
% end

% Combine for writing
uvw_pcui(:,:,:,1) = u_pcui;
uvw_pcui(:,:,:,2) = v_pcui;
uvw_pcui(:,:,:,3) = w_pcui;

% Write PCUI binary files depending on the number of processors
write_binary_file_pcui(working_folder, fname_rho_full_to_PCUI, params, rho_full_pcui);
write_binary_file_pcui(working_folder, fname_UVW_to_PCUI, params, uvw_pcui); 

%%
tmp = isnan(u_pcui(:));

figure;
hold on;
plot3(x_pcui(:),y_pcui(:),z_pcui(:),'k.');
plot3(x_pcui(tmp),y_pcui(tmp),z_pcui(tmp),'r.');

sum(isnan(u_pcui(:)))
sum(isnan(v_pcui(:)))
sum(isnan(rho_full_pcui(:)))



