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
isall = 0;
count = 1;
skip_vector = 8;


iskip =1;
istart =  1;
iend = 2400;
working_folder = '/work/knelson3/videoForOliver';

fname_xyz = 'xyz';
fname_uvw = 'output_UVW';

% -------------------------------------------------------------------------
% Get problem parameters and variables from PCUI
% -------------------------------------------------------------------------
% read the file containing the parameter definition
%ftext = fileread(fullfile(working_folder, 'io.f'));
ftext = fileread(fullfile(working_folder, 'pcuiRunParams.txt'));
params.dt = variable_value_pcui('dtime',ftext);
params.molecular_viscosity = variable_value_pcui('vis',ftext);
params.eddy_viscosity = variable_value_pcui('ak',ftext);
params.nsteps = variable_value_pcui('nstep',ftext);
params.nsave = variable_value_pcui('nsave',ftext);
params.dpdxSteady = variable_value_pcui('dpdxSteady',ftext);
params.waveMag = variable_value_pcui('waveMag',ftext);
params.Twave = variable_value_pcui('Twave',ftext);
params.waves = variable_value_pcui('waves',ftext);
params.rho_knot = 1000;%variable_value_pcui('rhoWater',ftext);

% read the file containing the grid size and processor definitions
ftext = fileread(fullfile(working_folder, 'size.inc'));
params.ni = variable_value_pcui('ni',ftext);
params.nj = variable_value_pcui('nj',ftext);
params.nk = variable_value_pcui('nk',ftext);
params.px = variable_value_pcui('px',ftext);
params.py = variable_value_pcui('py',ftext);
params.pz = variable_value_pcui('pz',ftext);

% read the files containing the grid definition and assemble into array
% (does not include ghost grid points)
[x,y,z] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,isall,2);

'i have x,y,z data'
params

% -------------------------------------------------------------------------
% Plot in 2D
% -------------------------------------------------------------------------
% Get the slice in the middle
middlez = floor(length(z(1,1,:))/2);
middlex = floor(length(z(:,1,1))/2);
middley = floor(length(z(1,:,1))/2);


x_xyplane(:,:) = squeeze(x(:,:,middlez));
y_xyplane(:,:) = squeeze(y(:,:,middlez));
z_xyplane(:,:) = squeeze(z(:,:,middlez));
clear x y z

count =1;
for istep = istart:iskip:iend
    istep 
    % Read velocity field
    [u,v,w] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
        params, 1,isall,2);


    u_xyplane(:,:,count) = squeeze(u(:,:,middlez));
    v_xyplane(:,:,count) = squeeze(v(:,:,middlez));
    w_xyplane(:,:,count) = squeeze(w(:,:,middlez));
    
  %  u_yzplane(:,:,count) = squeeze(u(middlex,:,:));
  %  v_yzplane(:,:,count) = squeeze(v(middlex,:,:));
  %  w_yzplane(:,:,count) = squeeze(w(middlex,:,:));
    
  %  u_xzplane(:,:,count) = squeeze(u(:,middley,:));
  %  v_xzplane(:,:,count) = squeeze(v(:,middley,:));
  %  w_xzplane(:,:,count) = squeeze(w(:,middley,:));
    
    velMag(:,:,count)=sqrt(u_xyplane(:,:,count).^2+w_xyplane(:,:,count).^2+v_xyplane(:,:,count).^2);
    count = count +1; 
end

save([working_folder '/velPlanes.mat'],'u_xyplane','v_xyplane','w_xyplane','velMag','x_xyplane','y_xyplane','z_xyplane','params')
