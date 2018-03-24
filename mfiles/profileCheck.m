close all; clear all
addpath('./Functions');
istart =1;
iskip = 1;

working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/strat500_1';
%working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/waveWithCurrents/dataWC350_19';
nu = 10^-6;

% These are the output files with the filenames stripped out of extensions
% (extensions are chosen automatically based on number of processors).
fname_drive = 'outputpVal_drive';
fname_umean = 'outputp_umean';
fname_vmean = 'outputp_vmean';
fname_wmean = 'outputp_wmean';
fname_uTurb = 'outputp_uTurb';
fname_vTurb = 'outputp_vTurb';
fname_wTurb = 'outputp_wTurb';
fname_uvRey = 'outputp_uvRey';
fname_uwRey = 'outputp_uwRey';
fname_vwRey = 'outputp_vwRey';
fname_phase = 'outputpVal_wavephase';


fname_Sed = 'outputp_cSed';
fname_vCsed = 'outputp_vCsed';
fname_sedTotal1 = 'outputpVal_sedTotal1';
fname_sedTotal2 = 'outputpVal_sedTotal2';

% -------------------------------------------------------------------------
% Get problem parameters and variables from PCUI
% -------------------------------------------------------------------------
% read the file containing the parameter definition
ftext = fileread(fullfile(working_folder, 'pcuiRunParams.txt'));
params.dt = variable_value_pcui('dtime',ftext);
params.molecular_viscosity = variable_value_pcui('vis',ftext);
params.eddy_viscosity = variable_value_pcui('ak',ftext);
params.nsteps = variable_value_pcui('nstep',ftext);
params.nsave = variable_value_pcui('nsave',ftext);
params.dpdxSteady = variable_value_pcui('dpdxSteady',ftext);
params.ieddy = variable_value_pcui('ieddy',ftext);
params.waves = variable_value_pcui('waves',ftext);
params.rho_knot = variable_value_pcui('rhoWater',ftext);
params.ised = variable_value_pcui('ised',ftext);
params.ws = variable_value_pcui('ws',ftext);

% read the file containing the grid size and processor definitions
ftext = fileread(fullfile(working_folder, 'size.inc'));
params.ni = variable_value_pcui('ni',ftext);
params.nj = variable_value_pcui('nj',ftext);
params.nk = variable_value_pcui('nk',ftext);
params.px = variable_value_pcui('px',ftext);
params.py = variable_value_pcui('py',ftext);
params.pz = variable_value_pcui('pz',ftext);
f = dir([working_folder '/' fname_wmean '.100']);
istop = f.bytes/(8.1250*params.nj/params.py)-5;
count = 1;
for istep = istart:iskip:istop
    umean(:,count) = read_binary_file_pcui_profile(working_folder, fname_umean, istep,params,0,0);
    vmean(:,count) = read_binary_file_pcui_profile(working_folder, fname_vmean, istep,params,0,0);
    wmean(:,count) = read_binary_file_pcui_profile(working_folder, fname_wmean, istep,params,0,0);
    uTurb(:,count) = read_binary_file_pcui_profile(working_folder, fname_uTurb, istep,params,0,0);
    vTurb(:,count) = read_binary_file_pcui_profile(working_folder, fname_vTurb, istep,params,0,0);
    wTurb(:,count) = read_binary_file_pcui_profile(working_folder, fname_wTurb, istep,params,0,0);
    uvRey(:,count) = read_binary_file_pcui_profile(working_folder, fname_uvRey, istep,params,0,0);
    uwRey(:,count) = read_binary_file_pcui_profile(working_folder, fname_uwRey, istep,params,0,0);
    vwRey(:,count) = read_binary_file_pcui_profile(working_folder, fname_vwRey, istep,params,0,0);
    dpdx(1,count) = read_binary_file_pcui_value(working_folder,   fname_drive, istep,1);
    phase(1,count) = read_binary_file_pcui_value(working_folder,  fname_phase, istep,1);
    if params.ised == 1
        Csed(:,count) = read_binary_file_pcui_profile(working_folder, fname_Sed, istep,params,0,0);
        Csed(:,count) = read_binary_file_pcui_profile(working_folder, fname_Sed, istep,params,0,0);
        CsedTotal1(1,count) = read_binary_file_pcui_value(working_folder,   fname_sedTotal1, istep,1);
        CsedTotal2(1,count) = read_binary_file_pcui_value(working_folder,   fname_sedTotal2, istep,1);
    end
    count =count+1;
end
%%%%%%%%%%%%%%%% Domain 2%%%%%%%%%%%%%%%%%%%%%%%%%%%
NJ=8*16;
NK=20*16;
NI=48*16;
px = 24;
py = 2;
pz = 10;
max_aspect_ratio=8.2;
num_wall_unit_bed = 1.2;
stretch_factor = 5/100;
u_star = .01; %friction velocity
wall_unit = nu/u_star; %length of wall unit
dy_min = wall_unit*num_wall_unit_bed; %dy at bed
dx= dy_min*max_aspect_ratio*2;
dz= dy_min*max_aspect_ratio;

stop_stretch = floor(log(dz/dy_min)/log(1+stretch_factor)+1); %cell number where vertical stretching stops
%dy of each cell
for i= 1:NJ
    if i==1
        dy(i)=dy_min;
    elseif i>stop_stretch
        dy(i)=dz;
    else
        dy(i)=dy_min*(1+stretch_factor)^(i-1);
    end
    aspect_ratio(i)=dx/dy(i);
end

%y coordinate of cell center (note: y is vertical coordinate in PCUI)
y(1)=dy_min/2;
for i= 1:NJ-1
    growth_rate(i) = dy(i+1)/dy(i);
    y(i+1)=y(i)+(dy(i)+dy(i+1))/2;
end
H = sum(dy);
u_star_plot = sqrt(params.dpdxSteady*H/params.rho_knot);

numFrames = length(CsedTotal1);

fig1 = figure;
set(0,'DefaultAxesFontSize',16)
set(fig1, 'Position', [100, 100, 1049, 895]);
ha = tight_subplot(1,2,[.05 .05],[.1 .04],[.1 .08]);
count =1;
for i = 1:numFrames
    %axes(ha(1))
    cla
    plot(umean(:,i)',y,'b',wmean(:,i)',y,'r',vmean(:,i)',y,'g','linewidth',2)
    hold on
    ylabel('vertical position (m)')
    xlabel('u_{mean} (m/s)')
    ylabel('vertical position (m)')
    %set(gca,'ylim',[0 .01])
    
    leg1 = legend('streamwise','cross-stream','vertical');
    set(leg1,'interpreter','Latex','FontSize',16,'FontWeight','bold')
    
    hold on
    drawnow;
    %pause
    %delete(timelabel)
end

if params.ised == 1
    fig01 = figure;
    set(0,'DefaultAxesFontSize',16)
    set(fig01, 'Position', [100, 100, 1049, 895]);
    ha = tight_subplot(1,2,[.05 .05],[.1 .04],[.1 .08]);
    
    count =1;
    for i = 1:numFrames
        axes(ha(1))
        cla
        plot(uTurb(:,i)'/u_star_plot,y,'b',wTurb(:,i)'/u_star_plot,y,'r',vTurb(:,i)'/u_star_plot,y,'g','linewidth',2)
        
        legend('streamwise','cross-stream','vertical','postion','northeast')
        hold on
        xlab1 = xlabel('$\frac{u_{i,rms}}{u_{*}}$');
        set(xlab1,'interpreter','Latex','FontSize',16,'FontWeight','bold')
        ylab = ylabel('vertical position (m)');
        set(ylab,'interpreter','Latex','FontSize',16,'FontWeight','bold')
        set(gca,'ylim',[0 H+0.05*H], 'xlim',[0 max(max(uTurb))/u_star_plot])
        
        axes(ha(2))
        cla
        plot(Csed(:,i)*1000,y,'k','linewidth',2)
        xlab1 = xlabel('$C \, (mg/L)$');
        set(xlab1,'interpreter','Latex','FontSize',16,'FontWeight','bold')
        %set(gca,'xlim',[0 1000],'ylim',[0 H+0.05*H],'yticklabel',[])
        set(gca,'xlim',[0 1000],'ylim',[0 H],'yticklabel',[])
        
        drawnow;
    end
end

umeanModel = umean;
yModel = y;
save('1stStep.mat','umeanModel','yModel');
