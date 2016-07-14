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

smooth = 1; % Interpolate shading in plot
isprint = 0; % Print plot to eps file
plot_name = 'channel flow';
moviename = 'blob_SCCandIso.gif';
isall = 0;
movie = true;
count = 1;
skip_vector = 4;
load('SedCmap.mat')

%%%%%% Testing Simulations
working_folder = '~/Desktop/SettlingVelocityTest/Circle_20000_fine';
%working_folder = '~/Desktop/PCUI_repository/pcui-3d/examples/sedModel';
%working_folder = 'F:\PCUI Results\SettlingVelocityTest\1000mgL';
%working_folder = 'C:\Users\knelson81\Desktop\PCUI_Data Share\Simulations\sedModel';


% These are the output files with the filenames stripped out of extensions
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
fname_rho = 'output_rho';
fname_uvw = 'output_UVW';
fname_Csed = 'output_Csed';
fname_pressure = 'output_p';

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
params.dpdxSteady = variable_value_pcui('dpdxSteady',ftext);
params.waveMag = variable_value_pcui('waveMag',ftext);
params.Twave = variable_value_pcui('Twave',ftext);
params.waves = variable_value_pcui('waves',ftext);

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



% Theoretical Solution for Steady Current
%y_theo = linspace(params.by/(params.nj*2),params.by-params.by/(params.nj*2),params.nj);
y_theo = linspace(0,params.by,params.nj);
u_theo = 1/(2*params.molecular_viscosity)*-params.dpdxSteady*y_theo.^2-1/(params.molecular_viscosity)*-params.dpdxSteady*params.by*y_theo;

% Variables for theoretical wave solution
omega = 2*pi/params.Twave;
delta = sqrt(omega/(2*params.molecular_viscosity));

% read the files containing the grid definition and assemble into array
% (does not include ghost grid points)
[x,y,z] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,isall,2);
dx =x(2,1,1)-x(1,1,1);
dy =y(1,2,1)-y(1,1,1);
dz =z(1,1,2)-z(1,1,1);
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

% if ~movie
%     fig1 = figure;
% %    set(fig1,'Position',[100 100 900 900]);
%     set(fig1,'Position',[457 131 902 692]);
% 
%     hold
%     %clf
%     % set(fig1,'Renderer','zbuffer');
%     % set(fig1,'Color','black');
% end


%%For MATLAB 2015 gif creator
fig1 = figure;
colormap(cmap);
    set(fig1,'Position',[457 131 902 692]);
%set(fig1,'Position',[100 100 900 900]);
%clf

iskip = 1;
istart =  1;
iend = floor(params.nsteps/params.nsave);
for istep = istart:iskip:iend
    
    %%%%% For older MATLAB versions
%     if movie
%         close all;
%         fig1 = figure;
%         set(fig1,'Position',[100 100 900 900]);
%         hold
%     end
    
    if(params.waves == 1)
        u_theo = params.waveMag*cos(omega*(istep-1)*params.dt*params.nsave)-params.waveMag*exp(-y_theo/delta).*cos(omega*(istep-1)*params.dt*params.nsave-y_theo/delta);
        max(u_theo)
    end
    
    % Read velocity field
    [u,v,w] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
        params, 1,isall,2);
    u_xyplane = squeeze(u(:,:,middlez));
    v_xyplane = squeeze(v(:,:,middlez));
    w_xyplane = squeeze(w(:,:,middlez));
    
    u_yzplane = squeeze(u(middlex,:,:));
    v_yzplane = squeeze(v(middlex,:,:));
    w_yzplane = squeeze(w(middlex,:,:));
    
    u_xzplane = squeeze(u(:,middley,:));
    v_xzplane = squeeze(v(:,middley,:));
    w_xzplane = squeeze(w(:,middley,:));
    
    % Read density field
    rho = read_binary_file_pcui(working_folder, fname_rho, istep, ...
        params, 0,isall,2);
    rho = squeeze(rho(:,:,middlez));
    
    % Read Sediment field
    Csed = read_binary_file_pcui(working_folder, fname_Csed, istep, ...
        params, 0,isall,2);
    MassSed(istep)=sum(sum(sum(Csed)))*dx*dy*dz;
    Csed_all = Csed;
    Csed = squeeze(Csed(:,:,middlez));
    
    % Read Sediment field
    pressure = read_binary_file_pcui(working_folder, fname_pressure, istep, ...
        params, 0,isall,1);
    pressure = squeeze(pressure(:,:,middlez));
    
%     % Plot
%     subplot(2,2,1)
%     cla
%     h1=pcolor(x_xyplane(2:end-1,2:end-1),y_xyplane(2:end-1,2:end-1),pressure);
%     %set(h1, 'EdgeColor', 'none');
%     %h1=quiver(x_xyplane(1:skip_vector:end,1:skip_vector:end),y_xyplane(1:skip_vector:end,1:skip_vector:end),u_xyplane(1:skip_vector:end,1:skip_vector:end),v_xyplane(1:skip_vector:end,1:skip_vector:end));
%     xlabel('x (m)')
%     ylabel('y (m)')
%     title('Pressure in xy-plane')
%     axis([0 params.bx 0 params.by])
%     if smooth
%         shading interp
%     end
%     colorbar
%     %caxis([-1 1])
%     annotation(gcf,'textbox',...
%     [0.421 0.934444444444444 0.026777777777778 0.0288888888888889],...
%     'String',{'Pa'},...
%     'LineStyle','none',...
%     'FitBoxToText','off');
%     
%     subplot(2,2,2)
%     cla
%     %plot(u_xyplane(middlex,:),y_xyplane(middlex,:),'k',w_xyplane(middlex,:),y_xyplane(middlex,:),'r',v_xyplane(middlex,:),y_xyplane(middlex,:),'b')
%     %plot(u_xyplane(middlex,:),y_xyplane(middlex,:),'k')
%     %plot(u_xyplane(middlex,:)-min(u_xyplane(middlex,:)),y_xyplane(middlex,:),'k')
%     h1=pcolor(x_xyplane,y_xyplane,sqrt(u_xyplane.^2+v_xyplane.^2+w_xyplane.^2));
%     xlabel('x (m)')
%     ylabel('y (m)')
%     title('velocity magnitude')
%     axis([0 params.bx 0 params.by])
%     if smooth
%         shading interp
%     end
%     colorbar
%     %xlabel('velocity (m/s)')
%     %ylabel('vertical position (m)')
%     %title('vertical profiles')
%     %legend('stream wise','cross stream','vertical','location', 'southeast')
%     %set(gca,'xlim',[0 5*10^-7])
%     %set(gca,'xlim',[-5*10^-7 5*10^-7])
%     % set(gca,'ylim',[0 params.by])
%     hold on
%     h1=quiver(x_xyplane(1:skip_vector:end,1:skip_vector:end),y_xyplane(1:skip_vector:end,1:skip_vector:end),...
%         u_xyplane(1:skip_vector:end,1:skip_vector:end),v_xyplane(1:skip_vector:end,1:skip_vector:end),'w');
%     annotation(gcf,'textbox',...
%     [0.862111111111111 0.941111111111111 0.0134444444444445 0.0244444444444445],...
%     'String','m/s',...
%     'LineStyle','none',...
%     'FitBoxToText','off');
%     
%     subplot(2,2,3)
%     cla
%     h1=pcolor(x_xyplane,y_xyplane,Csed*1000);
%     %set(h1, 'EdgeColor', 'none');
%     xlabel('x (m)')
%     ylabel('y (m)')
%     title('C_{sed} xy-plane')
%     %axis([0 params.bx 0 params.by])
%     if smooth
%         shading interp
%     end
%     colorbar
%     hold on
%     h1=quiver(x_xyplane(1:skip_vector:end,1:skip_vector:end),y_xyplane(1:skip_vector:end,1:skip_vector:end),...
%         u_xyplane(1:skip_vector:end,1:skip_vector:end),v_xyplane(1:skip_vector:end,1:skip_vector:end),'w');
%     %caxis([0 100])
%     annotation(gcf,'textbox',...
%     [0.431 0.472222222222222 0.026777777777778 0.0244444444444445],...
%     'String','mg/L',...
%     'LineStyle','none',...
%     'FitBoxToText','off');
%     
%     subplot(2,2,4)
%     cla
%     h1=pcolor(x_xyplane,y_xyplane,rho);
%     %set(h1, 'EdgeColor', 'none');
%     xlabel('x (m)')
%     ylabel('y (m)')
%     title('\rho xy-plane')
%     %axis([0 params.bx 0 params.by])
%     if smooth
%         shading interp
%     end
%     colorbar
%     annotation(gcf,'textbox',...
%     [0.886555555555556 0.471111111111111 0.0356666666666667 0.0244444444444445],...
%     'String','kg/m^3',...
%     'LineStyle','none',...
%     'FitBoxToText','off');
%     
%     %     subplot(2,2,4)
%     %     cla
%     %     plot(Csed(middlex,:),y_xyplane(middlex,:),'Color',[.7 .5 0])
%     %     xlabel('C_{sed} (mg/L)')
%     %     ylabel('vertical position (m)')
%     %     title('Vertical Sediment Concentration')
%     %     set(gca,'xlim',[0 1.2])
%     %     set(gca,'ylim',[0 params.by])
%     %
%     %
%     tim = (istep-1)*params.dt*params.nsave;
%     str = [num2str(tim) '/' num2str(params.dt*params.nsteps) ' s'];
%     timelabel = annotation(gcf,'textbox',[0.471 0.556 0.098 0.020000000000001],...
%         'String',str,...
%         'FitBoxToText','off',...
%         'LineStyle','none');

    subplot(1,2,1)
    cla
    h1=pcolor(x_xyplane,y_xyplane,Csed*1000);
    %set(h1, 'EdgeColor', 'none');
    xlabel('x (m)')
    ylabel('y (m)')
    title('C_{sed} xy-plane')
    %axis([0 params.bx 0 params.by])
    if smooth
        shading interp
    end
    axis equal
    colormap(cmap);
    colorbar
    hold on
    h1=quiver(x_xyplane(1:skip_vector:end,1:skip_vector:end),y_xyplane(1:skip_vector:end,1:skip_vector:end),...
        u_xyplane(1:skip_vector:end,1:skip_vector:end),v_xyplane(1:skip_vector:end,1:skip_vector:end),'w');
    %caxis([0 100])
    annotation(gcf,'textbox',...
    [0.480889135254989 0.535806037251124 0.0267777777777781 0.0244444444444445],...
    'String','mg/L',...
    'LineStyle','none',...
    'FitBoxToText','off');

    subplot(1,2,2)
    iso = patch(isosurface(Csed_all*1000,100));
    iso.FaceColor = [.7 .5 0];
    iso.EdgeColor = 'none';
    set(gca,'xlim',[1 128],'ylim',[1 64],'zlim',[1 64])
    view(gca,[89.7 -74]);
    camlight
    lighting gouraud
    axis equal
    set(gca,'xtick',[],'ytick',[],'ztick',[])
    drawnow;
    
    if movie;
        frame = getframe(gcf);
        %im = frame2im(frame);
        %[imind,cm] = rgb2ind(im,256);
        
        if count == 1
            [im,map] = rgb2ind(frame.cdata,256,'nodither');
            im(1,1,1,length(istart:iskip:iend)) = 0;
            %imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            im(:,:,1,count) = rgb2ind(frame.cdata,map,'nodither');
            %imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
        clf
        count=count+1;
    else
        drawnow;
        delete(timelabel)
    end
end

%
if movie
    imwrite(im,map,[working_folder '\' moviename],'DelayTime',.1,'LoopCount',inf);
end


% if isprint:q

%     file_name = fullfile(save_folder, plot_name);
%     export_fig(fig1,file_name, '-eps','-a1','-q101');
%     saveas(fig1,[file_name '.fig'],'fig');
% end