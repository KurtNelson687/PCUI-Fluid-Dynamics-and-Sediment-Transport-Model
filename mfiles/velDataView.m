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

'im in'

smooth = 1; % Interpolate shading in plot
isprint = 0; % Print plot to eps file
plot_name = 'channel flow';
loglaw = 1;
isall = 0;
count = 1;
skip_vector = 8;


iskip = 1;
istart =  1;
iend = 9;%floor(params.nsteps/params.nsave)

%%%%%% Testing Simulations

%working_folder = '/p/work1/knelson3/Wave4_2_smallresid';
working_folder = '/work/knelson3/Domain2';
%working_folder = '/work/knelson3/small_1';

fname_xyz = 'xyz';
fname_rho = 'output_rho';
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
% read the file containing the domain definition
ftext = fileread(fullfile(working_folder, 'cavity.f'));
params.wallUnit = variable_value_pcui('wallUnit',ftext);
params.maxRatio = variable_value_pcui('maxRatio',ftext);
params.numWallUnits = variable_value_pcui('numWallUnits',ftext);

params.bx = params.wallUnit*params.numWallUnits*params.maxRatio*params.ni;
params.bz = params.wallUnit*20*params.nk;
params.by = y(1,end,1)+params.wallUnit*params.numWallUnits*params.maxRatio/2;

params

% Theoretical Solution for Steady Current
y_theo = linspace(params.by/(params.nj*2),params.by-params.by/(params.nj*2),params.nj);
y_theo = linspace(0,params.by,params.nj);
u_theo = 1/(2*params.molecular_viscosity)*-params.dpdxSteady*y_theo.^2-1/(params.molecular_viscosity)*-params.dpdxSteady*params.by*y_theo;

% Variables for theoretical wave solution
omega = 2*pi/params.Twave;
delta = sqrt(omega/(2*params.molecular_viscosity));

% Log-law
u_star = sqrt(params.dpdxSteady*y(1,end,1)/params.rho_knot);
zo = params.molecular_viscosity/(9*u_star);
u_log = u_star/0.41*log(y_theo/zo);

% -------------------------------------------------------------------------
% Plot in 2D
% -------------------------------------------------------------------------
% Get the slice in the middle
middlez = floor(length(z(1,1,:))/2);
middlex = floor(length(z(:,1,1))/2);
middley = floor(length(z(1,:,1))/2);

x_xyplane = squeeze(x(:,:,middlez));
y_xyplane = squeeze(y(:,:,middlez));
z_xyplane = squeeze(z(:,:,middlez));

x_yzplane = squeeze(x(middlex,:,:));
y_yzplane = squeeze(y(middlex,:,:));
z_yzplane = squeeze(z(middlex,:,:));

x_xzplane = squeeze(x(:,middley,:));
y_xzplane = squeeze(y(:,middley,:));
z_xzplane = squeeze(z(:,middley,:));

clear x y z

count =1;
for istep = istart:iskip:iend
    if(params.waves == 1)
        u_theo = params.waveMag*cos(omega*(istep-1)*params.dt*params.nsave)-params.waveMag*exp(-y_theo/delta).*cos(omega*(istep-1)*params.dt*params.nsave-y_theo/delta);
        max(u_theo)
    end
    istep 
    % Read velocity field
    [u,v,w] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
        params, 1,isall,2);

    u_mean(count,:) = mean(mean(u,1),3);

    u_xyplane(:,:,count) = squeeze(u(:,:,middlez));
    v_xyplane(:,:,count) = squeeze(v(:,:,middlez));
    w_xyplane(:,:,count) = squeeze(w(:,:,middlez));
    
    u_yzplane(:,:,count) = squeeze(u(middlex,:,:));
    v_yzplane(:,:,count) = squeeze(v(middlex,:,:));
    w_yzplane(:,:,count) = squeeze(w(middlex,:,:));
    
    u_xzplane(:,:,count) = squeeze(u(:,middley,:));
    v_xzplane(:,:,count) = squeeze(v(:,middley,:));
    w_xzplane(:,:,count) = squeeze(w(:,middley,:));
    
    velMag(:,:,count)=sqrt(u_xyplane(:,:,count).^2+w_xyplane(:,:,count).^2+v_xyplane(:,:,count).^2);
    count = count +1; 
end

save([working_folder '/velPlanes.mat'],'u_xyplane','v_xyplane','w_xyplane','u_yzplane','v_yzplane','w_yzplane','u_xzplane',...
        'v_xzplane','w_xzplane','velMag','x_xyplane','y_xyplane','z_xyplane','x_xzplane','y_xzplane','z_xzplane','u_mean','params')
