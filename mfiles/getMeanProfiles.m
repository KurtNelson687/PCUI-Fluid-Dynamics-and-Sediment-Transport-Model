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
loglaw = 1;
moviename = 'turbulentStatis.gif';
isall = 0;
movie = true;
count = 1;
skip_vector = 10;

%%%%%% Testing Simulations
working_folder = '~/Desktop/PCUI_repository/pcui-3d/examples/sedModel';

%working_folder = 'C:\Users\knelson81\Desktop\PCUI_Data Share\Simulations\sedModel';
% working_folder = 'C:\Users\knelson81\Desktop\PCUI_Data Share\Simulations\currentsAndWaves';
% working_folder = 'C:\Users\knelson81\Desktop\PCUI_Data Share\Simulations\couetteFlow';
% working_folder = 'C:\Users\knelson81\Desktop\PCUI_Data Share\Simulations\PoiseuilleFlow';

%working_folder = 'C:\Users\knelson81\Desktop\PCUI_Data Share\Simulations\wavesOnly';
% working_folder = 'C:\Users\knelson81\Desktop\PCUI_Data Share\Simulations\wavesOnlyWithLidBC';
% working_folder = 'C:\Users\knelson81\Desktop\PCUI_Data Share\Simulations\wavesOnlyWithPressure';
% working_folder = 'C:\Users\knelson81\Desktop\PCUI_Data Share\Simulations\sedimentModel';


%save_folder =    'C:\Users\knelson81\Desktop\PCUI_Data Share\Simulations\pcui_steady_current';


% %%%% Current only
% working_folder = 'C:\Users\knelson81\Dropbox\Data\PCUI Transfer\Simulations\Viscosity_1\currentOnly';
% save_folder =    'C:\Users\knelson81\Dropbox\Data\PCUI Transfer\Simulations\Viscosity_1\currentOnly';

%%%% Waves only
% working_folder = 'C:\Users\knelson81\Dropbox\Data\PCUI Transfer\Simulations\Viscosity_1\wavesOnly';
% save_folder =    'C:\Users\knelson81\Dropbox\Data\PCUI Transfer\Simulations\Viscosity_1\wavesOnly';

%%%% Waves and currents
% working_folder = 'C:\Users\knelson81\Dropbox\Data\PCUI Transfer\Simulations\Viscosity_1\wavesAndCurrents';
% save_folder =    'C:\Users\knelson81\Dropbox\Data\PCUI Transfer\Simulations\Viscosity_1\wavesAndCurrents';

% These are the output files with the filenames stripped out of extensions
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
%fname_rho = 'output_S';
fname_rho = 'output_rho';

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

% read the file containing the domain definition
ftext = fileread(fullfile(working_folder, 'cavity.f'));
params.wallUnit = variable_value_pcui('wallUnit',ftext);
params.maxRatio = variable_value_pcui('maxRatio',ftext);
params.numWallUnits = variable_value_pcui('numWallUnits',ftext);

params.bx = params.wallUnit*params.numWallUnits*params.maxRatio*params.ni;
params.bz = params.wallUnit*params.numWallUnits*params.maxRatio*params.nk;
params.by = y(1,end,1)+params.wallUnit*params.numWallUnits*params.maxRatio/2;


% Theoretical Solution for Steady Current
%y_theo = linspace(params.by/(params.nj*2),params.by-params.by/(params.nj*2),params.nj);
y_theo = y(1,:,1);%linspace(0,params.by,params.nj);
u_theo = 1/(2*params.molecular_viscosity)*-params.dpdxSteady*y_theo.^2-1/(params.molecular_viscosity)*-params.dpdxSteady*params.by*y_theo;

% Variables for theoretical wave solution
omega = 2*pi/params.Twave;
delta = sqrt(omega/(2*params.molecular_viscosity));

% Log-law
u_star = sqrt(params.dpdxSteady*(y(1,end,1)+(y(1,end,1)-y(1,end-1,1))/2)/params.rho_knot);
zo = params.molecular_viscosity/(9*u_star);
u_log = u_star/0.41*log(y_theo/zo);
% figure
% plot(u_theo,y_theo,'k',u_log,y_theo,'b')
% xlabel('velocity (m/s)')
% ylabel('depth (m)')
% legend('laminar','log law')

% -------------------------------------------------------------------------
% Plot in 2D
% -------------------------------------------------------------------------
% Get the slice in the middle
bottomy = find(y(1,:,1)>params.wallUnit*15,1,'first');

middlez = find(z(1,1,:)>params.bz/2,1,'first');
middlex = find(x(:,1,1)>params.bx/2,1,'first');
middley = find(y(1,:,1)>params.by/2,1,'first');

topy = find(y(1,:,1)>params.by*3/4,1,'first');



x_xyplane = squeeze(x(:,:,middlez));
y_xyplane = squeeze(y(:,:,middlez));
z_xyplane = squeeze(z(:,:,middlez));

% x_yzplane = squeeze(x(middlex,:,:));
% y_yzplane = squeeze(y(middlex,:,:));
% z_yzplane = squeeze(z(middlex,:,:));
% 
% x_xzplane = squeeze(x(:,middley,:));
% y_xzplane = squeeze(y(:,middley,:));
% z_xzplane = squeeze(z(:,middley,:));

fig1 = figure(1);
clf
% set(fig1,'Renderer','zbuffer');
% set(fig1,'Color','black');
set(fig1,'Position',[200 200 1000 1000]);

iskip = 1;
istart =  1;
iend = 360;%floor(params.nsteps/params.nsave);

if ~movie
    fig1 = figure;
    set(fig1,'Position',[100 100 900 900]);
    hold
    %clf
    % set(fig1,'Renderer','zbuffer');
    % set(fig1,'Color','black');
end

for istep = istart:iskip:iend
    if(params.waves == 1)
        u_theo = params.waveMag*cos(omega*(istep-1)*params.dt*params.nsave)-params.waveMag*exp(-y_theo/delta).*cos(omega*(istep-1)*params.dt*params.nsave-y_theo/delta);
        max(u_theo)
    end
    
    % Read velocity field
    [u,v,w] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
        params, 1,isall,2);
    
    [m,n,p] = size(u);
    
    u_xyplane = squeeze(u(:,:,middlez));
     v_xyplane = squeeze(v(:,:,middlez));
%     w_xyplane = squeeze(w(:,:,middlez));
%     
%     u_yzplane = squeeze(u(middlex,:,:));
%     v_yzplane = squeeze(v(middlex,:,:));
%     w_yzplane = squeeze(w(middlex,:,:));
%     
%     u_xzplane = squeeze(u(:,middley,:));
%     v_xzplane = squeeze(v(:,middley,:));
%     w_xzplane = squeeze(w(:,middley,:));
    
    u_mean = mean(mean(u,1),3);
    v_mean = mean(mean(v,1),3);
    w_mean = mean(mean(w,1),3);
    
    u_meanMatrix = repmat(u_mean,[m,1,p]);
    w_meanMatrix = repmat(w_mean,[m,1,p]);
    v_meanMatrix = repmat(v_mean,[m,1,p]);
    
    u_prime = u-u_meanMatrix;
    v_prime = v-v_meanMatrix;
    w_prime = w-w_meanMatrix;
    
    %Turbulent Intensity
    u_turbIntenseity = (u_prime.^2).^0.5;
    v_turbIntenseity = (v_prime.^2).^0.5;
    w_turbIntenseity = (w_prime.^2).^0.5;
    u_turbIntenseity_mean(count,:) = mean(mean(u_turbIntenseity,1),3);
    v_turbIntenseity_mean(count,:) = mean(mean(v_turbIntenseity,1),3);
    w_turbIntenseity_mean(count,:) = mean(mean(w_turbIntenseity,1),3);
    
    %Reynolds Stresses
    uv = -1*u_prime.*v_prime;
    uv_mean(count,:) = mean(mean(uv,1),3);
    uw = -1*u_prime.*w_prime;
    uw_mean(count,:) = mean(mean(uw,1),3);  
    vw = -1*v_prime.*w_prime;
    vw_mean(count,:) = mean(mean(vw,1),3);  
        
    %Two point correlation functions (see Moin and Kim 1982)
    r=1:m;
    r3=1:p;
    R11_all_bottom(:,count) = TwoPointCorrelation(u_prime,r,bottomy,13,1); 
    R22_all_bottom(:,count) = TwoPointCorrelation(v_prime,r,bottomy,13,1); 
    R33_all_bottom(:,count) = TwoPointCorrelation(w_prime,r,bottomy,13,1);
    
    R11_all_cross_bottom(:,count) = TwoPointCorrelation(u_prime,r3,bottomy,13,3); 
    R22_all_cross_bottom(:,count) = TwoPointCorrelation(v_prime,r3,bottomy,13,3); 
    R33_all_cross_bottom(:,count) = TwoPointCorrelation(w_prime,r3,bottomy,13,3);
    
    R11_all_middle(:,count) = TwoPointCorrelation(u_prime,r,middley,13,1);
    R22_all_middle(:,count) = TwoPointCorrelation(v_prime,r,middley,13,1); 
    R33_all_middle(:,count) = TwoPointCorrelation(w_prime,r,middley,13,1);
    
    R11_all_top(:,count) = TwoPointCorrelation(u_prime,r,topy,13,1);
    R22_all_top(:,count) = TwoPointCorrelation(v_prime,r,topy,13,1);
    R33_all_top(:,count) = TwoPointCorrelation(w_prime,r,topy,13,1);
    
%     % Read density field
%     rho = read_binary_file_pcui(working_folder, fname_rho, istep, ...
%         params, 0,isall,2);
%     rho = squeeze(rho(:,:,middlez));
    
    % Plot
    subplot(2,2,1)
    cla
    h1=pcolor(x_xyplane,y_xyplane,u_xyplane);
    %set(h1, 'EdgeColor', 'none');
    xlabel('x (m)')
    ylabel('y (m)')
    title('u in xy-plane')
    axis([0 params.bx 0 params.by])
    if smooth
        shading interp
    end
    colorbar
    hold on
    h1=quiver(x_xyplane(1:skip_vector:end,1:skip_vector:end),y_xyplane(1:skip_vector:end,1:skip_vector:end),...
        u_xyplane(1:skip_vector:end,1:skip_vector:end),v_xyplane(1:skip_vector:end,1:skip_vector:end),'w');
    
    subplot(2,2,2)
    cla
    plot(u_log,y_theo,'--k',u_mean,y(1,:,1),'k',w_mean,y(1,:,1),'r',v_mean,y(1,:,1),'b','linewidth',2)
    xlabel('velocity (m/s)')
    ylabel('vertical position (m)')
    title('velocity profiles (m/s)')
    set(gca,'ylim',[0 params.by])
    hold on
    %plot(u(1,:,1),y(1,:,1), 'color', [0.5 0.5 0.5],'linewidth',0.5)
    legend('Theoretical Log Law','<streamwise>','<cross-stream>','<vertical>','location','northwest')
%     for i = 1:m
%         for k =1:p
%             plot(u(i,:,k),y(i,:,k), 'color', [0.5 0.5 0.5],'linewidth',0.5)
%         end
%     end
%     plot(u_log,y_theo,'--k',u_mean,y(1,:,1),'k','linewidth',3)
    
    subplot(2,2,3)
    cla
    plot(u_turbIntenseity_mean(count,:),y(1,:,1),'k',w_turbIntenseity_mean(count,:),y(1,:,1),'r',v_turbIntenseity_mean(count,:),y(1,:,1),'b','linewidth',2)
    xlabel('turbulent intensity (m/s)')
    ylabel('vertical position (m)')
    tit = title('$$<\sqrt{u_{i}^{''2}}>$$');
    set(tit,'Interpreter','latex')
    %set(gca,'ylim',[0 params.by])
    %hold on
    %plot(v(1,:,1),y(1,:,1), 'color', [0.5 0.5 0.5],'linewidth',0.5)
    legend('streamwise','cross-stream','vertical')
%     for i = 1:m
%         for k =1:p
%             plot(v(i,:,k),y(i,:,k), 'color', [0.5 0.5 0.5],'linewidth',0.5)
%         end
%     end
    %plot(v_mean,y(1,:,1),'k','linewidth',3


    subplot(2,2,4)
    cla
    plot(uv_mean(count,:),y(1,:,1),'m',uw_mean(count,:),y(1,:,1),'g',vw_mean(count,:),y(1,:,1),'y','linewidth',2)
    xlabel('(m^2/s^2)')
    ylabel('vertical position (m)')
    title('reynolds stresses')
    set(gca,'ylim',[0 params.by])
%     hold on
%     plot(w(1,:,1),y(1,:,1), 'color', [0.5 0.5 0.5],'linewidth',0.5)
    legend('-<u''v''>','-<u''w''>','-<v''w''>')
%     for i = 1:m
%         for k =1:p
%             plot(w(i,:,k),y(i,:,k), 'color', [0.5 0.5 0.5],'linewidth',0.5)
%         end
%     end
%     plot(w_mean,y(1,:,1),'k','linewidth',3)
   % set(legend1,...
    %    'Position',[0.43368817204301 0.480225127087872 0.135783945986497 0.0660856935366739]);
    
    tim = (istep-1)*params.dt*params.nsave;
    str = [num2str(tim) '/' num2str(params.dt*params.nsteps+180) ' s'];
    timelabel = annotation(gcf,'textbox',[0.471 0.556 0.098 0.020000000000001],...
        'String',str,...
        'FitBoxToText','off',...
        'LineStyle','none');
    
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
    else
        drawnow;
        delete(timelabel)
    end
    count=count+1;
end

if movie
    imwrite(im,map,[working_folder '/' moviename],'DelayTime',.2,'LoopCount',inf);
end

fig3 = figure;
set(0,'DefaultAxesFontSize',12)
fig3.PaperUnits = 'inches';
fig3.PaperPosition = [0 0 8.5 11];
fig3.PaperPositionMode = 'manual';
ha = tight_subplot(3,1,[.03 .03],[.1 .04],[.1 .08]);

sum_start=istep-100;
axes(ha(1))
R11_top = mean(R11_all_top(:,sum_start:count-1),2); 
R22_top = mean(R22_all_top(:,sum_start:count-1),2);
R33_top = mean(R33_all_top(:,sum_start:count-1),2);
plot(x(:,1,1)-x(1,1,1),R11_top,'*k',x(:,1,1)-x(1,1,1),R33_top,'or',x(:,1,1)-x(1,1,1),R22_top,'^b');
set(gca,'XTick', [],'Xlim',[0 params.bx]);

axes(ha(2))
R11_middle = mean(R11_all_middle(:,sum_start:count-1),2); 
R22_middle = mean(R22_all_middle(:,sum_start:count-1),2);
R33_middle = mean(R33_all_middle(:,sum_start:count-1),2);
plot(x(:,1,1)-x(1,1,1),R11_middle,'*k',x(:,1,1)-x(1,1,1),R33_middle,'or',x(:,1,1)-x(1,1,1),R22_middle,'^b');
set(gca,'XTick', [],'Xlim',[0 params.bx]);
ylabel('Two-point Correlation')

axes(ha(3))
R11_bottom = mean(R11_all_bottom(:,sum_start:count-1),2); 
R22_bottom = mean(R22_all_bottom(:,sum_start:count-1),2);
R33_bottom = mean(R33_all_bottom(:,sum_start:count-1),2);
plot(x(:,1,1)-x(1,1,1),R11_bottom,'*k',x(:,1,1)-x(1,1,1),R33_bottom,'or',x(:,1,1)-x(1,1,1),R22_bottom,'^b');
set(gca,'Xlim',[0 params.bx]);
xlabel('r_{1} (m)')
leg1 = legend('R^1_{11}','R^1_{33}','R^1_{22}');
set(leg1,...
    'Position',[0.82677448376384 0.883116885081871 0.043126684250536 0.0725940706290843]);
annotation(fig3,'textbox',...
    [0.101269541778977 0.124875124875125 0.140239892183289 0.0119880119880119],...
    'String','15 wall units from bed',...
    'LineStyle','none',...
    'FontSize',16,...
    'FitBoxToText','off');

annotation(fig3,'textbox',...
    [0.100730458221026 0.41958041958042 0.140239892183289 0.0119880119880119],...
    'String','1/2 channel depth',...
    'LineStyle','none',...
    'FontSize',16,...
    'FitBoxToText','off');

annotation(fig3,'textbox',...
    [0.101808625336928 0.716283716283716 0.140239892183289 0.0119880119880119],...
    'String',{'3/4 channel depth'},...
    'LineStyle','none',...
    'FontSize',16,...
    'FitBoxToText','off');

fig4 = figure;
set(0,'DefaultAxesFontSize',12)
plot(mean(u_turbIntenseity_mean(sum_start:count-1,:),1),y(1,:,1),'*k',mean(w_turbIntenseity_mean(sum_start:count-1,:),1),y(1,:,1),'or',mean(v_turbIntenseity_mean(sum_start:count-1,:),1),y(1,:,1),'^b')
ylabel('Vertical Position (m)')
xlabel('Turbulent Intensities (m/s)')
leg = legend('$<\sqrt{(u'')^2}>$','$<\sqrt{(w'')^2}>$','$<\sqrt{(v'')^2}>$');
set(leg,'interpreter','Latex','FontSize',12)

fig5 = figure;
set(0,'DefaultAxesFontSize',12)
plot(mean(uv_mean(sum_start:count-1,:),1),y(1,:,1),'*k',mean(uw_mean(sum_start:count-1,:),1),y(1,:,1),'or',mean(vw_mean(sum_start:count-1,:),1),y(1,:,1),'^b')
ylabel('Vertical Position (m)')
xlabel('Turbulent Intensities (m/s)')
leg = legend('$<-u''v''>$','$<-u''w''>$','$<-v''w''>$');
set(leg,'interpreter','Latex','FontSize',12)

fig6 = figure;
set(0,'DefaultAxesFontSize',12)
fig6.PaperUnits = 'inches';
fig6.PaperPosition = [0 0 8.5 11];
fig6.PaperPositionMode = 'manual';
ha = tight_subplot(2,1,[.03 .03],[.1 .04],[.1 .08]);

axes(ha(1))
R11_bottom_cross = mean(R11_all_cross_bottom(:,sum_start:count-1),2); 
R22_bottom_cross = mean(R22_all_cross_bottom(:,sum_start:count-1),2);
R33_bottom_cross = mean(R33_all_cross_bottom(:,sum_start:count-1),2);
plot(squeeze(z(1,1,:)-z(1,1,1)),R11_bottom_cross,'*k',squeeze(z(1,1,:)-z(1,1,1)),R33_bottom_cross,'or',squeeze(z(1,1,:)-z(1,1,1)),R22_bottom_cross,'^b');
set(gca,'XTick', [],'Xlim',[0 params.bz]);
leg2 = legend('R^2_{11}','R^2_{33}','R^2_{22}');



axes(ha(2))
plot(x(:,1,1)-x(1,1,1),R11_bottom,'*k',x(:,1,1)-x(1,1,1),R33_bottom,'or',x(:,1,1)-x(1,1,1),R22_bottom,'^b');
set(gca,'Xlim',[0 params.bx]);
xlabel('r_{1} (m)')
leg1 = legend('R^1_{11}','R^1_{33}','R^1_{22}');
set(leg1,...
    'Position',[0.82677448376384 0.883116885081871 0.043126684250536 0.0725940706290843]);
annotation(fig3,'textbox',...
    [0.101269541778977 0.124875124875125 0.140239892183289 0.0119880119880119],...
    'String','15 wall units from bed',...
    'LineStyle','none',...
    'FontSize',16,...
    'FitBoxToText','off');

annotation(fig3,'textbox',...
    [0.100730458221026 0.41958041958042 0.140239892183289 0.0119880119880119],...
    'String','1/2 channel depth',...
    'LineStyle','none',...
    'FontSize',16,...
    'FitBoxToText','off');

annotation(fig3,'textbox',...
    [0.101808625336928 0.716283716283716 0.140239892183289 0.0119880119880119],...
    'String',{'3/4 channel depth'},...
    'LineStyle','none',...
    'FontSize',16,...
    'FitBoxToText','off');



%     file_name = fullfile(save_folder, plot_name);
%     export_fig(fig1,file_name, '-eps','-a1','-q101');
%     saveas(fig1,[file_name '.fig'],'fig');
% end
