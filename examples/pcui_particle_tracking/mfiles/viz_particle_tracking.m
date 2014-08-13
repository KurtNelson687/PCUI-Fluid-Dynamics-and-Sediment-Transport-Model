% Displays PCUI particle tracking data 
clear all; clc; 

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
filename_XYZ = 'output_XYZ.dat';

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

% %read the file containing particle definitions
% ftext = fileread(fullfile(working_folder, 'particles.inc'));
% params.lbi = variable_value_pcui('lbi',ftext);
% params.ubi = variable_value_pcui('ubi',ftext);
% params.lbj = variable_value_pcui('lbj',ftext);
% params.ubj = variable_value_pcui('ubj',ftext);
% params.lbk = variable_value_pcui('lbk',ftext);
% params.ubk = variable_value_pcui('ubk',ftext);

%Load grid
[x,y,z] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,0);
x = squeeze(x(:,:,1));
y = squeeze(y(:,:,1));

%Load particle grid
[xp,yp,zp] = read_binary_particle_grid_pcui(working_folder, filename_XYZ, 1, params);
% xp = squeeze(xp(:,:,1));
% yp = squeeze(yp(:,:,1));
% zp = squeeze(yp(:,:,1));

%Find correct istep value
n = [0, params.nsave:params.nsave:params.nsteps, params.nsteps+1];

TEND = 600;
np = (40-17+1)*params.nj;
xpartall = nan(np,3,1);
for timestep = 0:params.nsave:TEND
    display(timestep);
    istep = find(n==timestep,1);
    xpartall(:,:,istep) = read_binary_particles_pcui(working_folder, filename_xpart, istep, np);
end

%%
clc; close all;
load '1b_bottom.mat';
% v = [34 44 101];
% v = [1 65 129 193];
v = [1 27 34 44 53 65 101 126 128 129 193 897 965 979];
figure;
set(gcf,'position',[400 400 1000 1000]);
hold on;
for timestep=0:params.nsave:TEND
    istep = find(n==timestep,1);
    display(istep);
    rho = read_binary_file_pcui(working_folder, fname_rho, istep, ...
                                 params, 0,0);  
%     [u,v,w] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
%                                  params, 1,0); 
    cla;
%     pcolor(x,y,squeeze(u(:,:,1))); shading flat; colorbar;
    contour(x,y,squeeze(rho(:,:,1)),[1 1],'r');
    plot(squeeze(xpartall(:,1,istep)),squeeze(xpartall(:,2,istep)),'k.');
%     plot(squeeze(xpartall(v,1,istep)),squeeze(xpartall(v,2,istep)),'r.');

    plot([0 0],[0 -params.by],'k-');
    plot([params.bx params.bx],[0 -0.14],'k-');
    plot([0 params.bx],[0 0],'k-');
    plot([0 1.675],[-params.by -params.by],'k-');
    plot(xbottom+1.675,bottom,'k-');
    axis equal;
    axis([-0.1 3.6 -0.6 0.005]); 
%     axis([2.5 3 -0.56 -0.15]);
%     plot3(xp(:),zp(:),yp(:),'b.');
%     plot3(squeeze(xpartall(:,1,istep)),squeeze(xpartall(:,3,istep)),squeeze(xpartall(:,2,istep)),'k.');
    drawnow;
    pause;
end
hold off;

%%
xgrid = read_binary_particles_pcui(working_folder, 'output_XYZ.dat', 0, params);
plot(squeeze(xgrid(:,1)),squeeze(xgrid(:,2)),'k.');
axis equal;