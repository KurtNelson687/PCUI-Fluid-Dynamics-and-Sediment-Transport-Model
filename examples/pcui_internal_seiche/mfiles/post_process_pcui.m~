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

smooth = 0; % Interpolate shading in plot
isprint = 0; % Print plot to eps file
plot_name = 'internal_seiche';
working_folder = '../';
save_folder =    '../figs';
% These are the output files with the filenames stripped out of extensions 
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
fname_rho = 'output_S';
fname_uvw = 'output_UVW';

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
% (does not include ghost grid points)
[x,y,z] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,0);

% -------------------------------------------------------------------------
% Plot in 2D
% -------------------------------------------------------------------------
% Get the slice in the middle
x = squeeze(x(:,:,floor(length(z(1,1,:)/2))));
y = squeeze(y(:,:,floor(length(z(1,1,:)/2))));
z = squeeze(z(:,:,floor(length(z(1,1,:)/2))));

fig1 = figure(1);
clf
set(fig1,'Renderer','zbuffer');
set(fig1,'Color','black');
set(fig1,'Position',[200 200 1000 1000]);

iskip = 2;
istart =  1;
iend = floor(params.nsteps/params.nsave);
for istep = istart:iskip:iend
    
    % Read velocity field
    [u,v,w] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
                                      params, 1,0);     
    u = squeeze(u(:,:,floor(length(z(1,1,:)/2))));
    v = squeeze(v(:,:,floor(length(z(1,1,:)/2))));
    w = squeeze(w(:,:,floor(length(z(1,1,:)/2))));
    
    % Read density field
    rho = read_binary_file_pcui(working_folder, fname_rho, istep, ...
                                 params, 0,0);     
    rho = squeeze(rho(:,:,floor(length(z(1,1,:)/2))));
    
    % Plot
    subplot(2,2,1)    
    cla    
    h1=pcolor(x,y,u);     
    %set(h1, 'EdgeColor', 'none');    
    xlabel('x (m)','color','w')
    ylabel('y (m)','color','w')
    title('u','color','w')
    axis([0 params.bx 0 params.by])
    if smooth
        shading interp
    end
    
    subplot(2,2,2)    
    cla    
    h2=pcolor(x,y,v);     
    %set(h2, 'EdgeColor', 'none');
    xlabel('x (m)','color','w')
    ylabel('y (m)','color','w')
    title('v','color','w')
    axis([0 params.bx 0 params.by])
    if smooth
        shading interp
    end
    
    subplot(2,2,3)    
    cla
    h3=pcolor(x,y,rho);     
    %set(h3, 'EdgeColor', 'none');   
    xlabel('x (m)','color','w')
    ylabel('y (m)','color','w')
    title('\rho','color','w')
    axis([0 params.bx 0 params.by])
    if smooth
        shading interp
    end
    %colorbar
    
    subplot(2,2,4)    
    cla    
    skipx = 4; intx = 1:skipx:params.ni;
    skipy = 4; inty = 1:skipy:params.nj;    
    h4=quiver(x(intx,inty),y(intx,inty),u(intx,inty),v(intx,inty));    
    xlabel('x (m)','color','w')
    ylabel('y (m)','color','w')
    title('u,v quiver plot','color','w')
    axis([0 params.bx 0 params.by])
    
    tim = (istep-1)*params.dt*params.nsave;
    str = [num2str(tim) '/' num2str(params.dt*params.nsteps) ' s'];
    uicontrol('Style', 'text', 'String', str, 'Units','normalized', ...
              'Position', [0.4 0.48 0.23 0.04]);
          
    drawnow;
end

if isprint
    file_name = fullfile(save_folder, plot_name);
    export_fig(fig1,file_name, '-eps','-a1','-q101');
    saveas(fig1,[file_name '.fig'],'fig');
end