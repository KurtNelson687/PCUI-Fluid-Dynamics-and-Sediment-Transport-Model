% Script: dataPlotting_velAndEquil.m
%
% Author: Kurt Nelson
%
% Purpose: This script creates 5 plots related to mean velocity profiles
% and checking for statistical equlibrium.The script requires processed 
% PCUI data from "dataExtractor.m".
%%
close all; clear all; clear all
addpath('./Functions');

%% Flags
savePlot = true; %save plots or not

%% Variables
data_folder = '/Users/kurtnelson/Desktop/Writing/WaveAndCurrent_sedFocus/mfiles/Data'; %PCUI data folder
data_fileWC = 'dataWCfineSpinup'; %PCUI processed data file for currents only
data_fileC = 'dataCurrentOnly.mat'; %PCUI processed data file for currents only
runsCup = [2,4,6]; % Runs to plot from wave and current simulaitons for stratified runs
runsCurrent = 1:3; % Runs to plot from current only simulaitons
runsAll = 1:7; % Runs to plot for wave and current simulations for stratified and unstratified runs
aveStart = 1; % Index to start phase averaging at
num_steps_period = 12*5; % Number of times data is recorded per wave cycle
phases = (0:num_steps_period)/num_steps_period*2*pi; % Number of phases to plot over
getDataAndGridInfo % Load PCUI data, plotting variables, grid data, and Moser and del Alamo data
plotSpecs.lineType = {'-','--',':'}; % Line types specific to this script
delta = sqrt(2*params.molecular_viscosity/(2*pi/3)); % Stokes layer thickness
kappa = 0.41; % von Karman constant
buffTop = 30; % z+ value chosen for top of the buffer layer
g = 9.81; % gravity m2/s

% x labels
xticklabels = {'0','\pi/4', '\pi/2','3\pi/4', '\pi','5\pi/4', '3\pi/2','7\pi/4', '2 \pi'};
xtick = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4,2*pi];

%%%%%%%%Compute Shear and density gradients%%%%%%%%%%

for plotNum = runsAll
    [m,n] = size(dataWC(plotNum).umean);
    c = zeros(m,5);
    for j = 3:m-2
        c(j,:) = fdcoeffF(1,y(j),y(j-2:j+2));
    end
    
    for j=1:m
        for time=1:n
            if j>2 && j<m-1
                dataWC(plotNum).Ushear(j,time) = c(j,:)*dataWC(plotNum).umean(j-2:j+2,time);
                dataWC(plotNum).Vshear(j,time) = c(j,:)*dataWC(plotNum).vmean(j-2:j+2,time);
                dataWC(plotNum).Wshear(j,time) = c(j,:)*dataWC(plotNum).wmean(j-2:j+2,time);
                dataWC(plotNum).rhoGrad(j,time) = c(j,:)*dataWC(plotNum).rhoMean(j-2:j+2,time);
            elseif j<=2
                dataWC(plotNum).Ushear(j,time) = (dataWC(plotNum).umean(j+1,time)-dataWC(plotNum).umean(j,time))/(y(j+1)-y(j));
                dataWC(plotNum).Vshear(j,time) = (dataWC(plotNum).vmean(j+1,time)-dataWC(plotNum).vmean(j,time))/(y(j+1)-y(j));
                dataWC(plotNum).Wshear(j,time) = (dataWC(plotNum).wmean(j+1,time)-dataWC(plotNum).wmean(j,time))/(y(j+1)-y(j));
                dataWC(plotNum).rhoGrad(j,time) = (dataWC(plotNum).rhoMean(j+1,time)-dataWC(plotNum).rhoMean(j,time))/(y(j+1)-y(j));
            else
                dataWC(plotNum).Ushear(j,time) = (dataWC(plotNum).umean(j,time)-dataWC(plotNum).umean(j-1,time))/(y(j)-y(j-1));
                dataWC(plotNum).Vshear(j,time) = (dataWC(plotNum).vmean(j,time)-dataWC(plotNum).vmean(j-1,time))/(y(j)-y(j-1));
                dataWC(plotNum).Wshear(j,time) = (dataWC(plotNum).wmean(j,time)-dataWC(plotNum).wmean(j-1,time))/(y(j)-y(j-1));
                dataWC(plotNum).rhoGrad(j,time) = (dataWC(plotNum).rhoMean(j,time)-dataWC(plotNum).rhoMean(j-1,time))/(y(j)-y(j-1));
            end
        end
    end
    
    for time=1:n
        %Production at cell center
        dataWC(plotNum).Production(:,time) = dataWC(plotNum).Ushear(:,time).*dataWC(plotNum).uvRey(:,time)+...
            dataWC(plotNum).Wshear(:,time).*dataWC(plotNum).vwRey(:,time)+dataWC(plotNum).Vshear(:,time).*dataWC(plotNum).vTurb(:,time).^2;
        
        %Bouyancy flux at cell center
        dataWC(plotNum).BuoyancyFlux(:,time) = g/params.rhoWater.*dataWC(plotNum).vCsed(:,time)*...
            (1-params.rhoWater/dataWC(plotNum).rhoSed);
        
        %Brunt Väisälä frequency at cell center
        dataWC(plotNum).BruntNCenter(:,time) = sqrt(-g/params.rhoWater*dataWC(plotNum).rhoGrad(:,time));
        
        %Gradient Richardson Number
        dataWC(plotNum).Ri(:,time) = dataWC(plotNum).BruntNCenter(:,time).^2./...
            (dataWC(plotNum).Ushear(j,time).^2+dataWC(plotNum).Wshear(j,time).^2);
    end
end

fig1 = figure;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = plotSpecs.fullWidthFigSize;
set(gca,'box','on')

yplus = [0.1:0.1:11.4, 12:10:500];
for i = 1:length(yplus)
    if yplus(i)<11.4
    theory(i) = yplus(i);
    else
    theory(i) = 1/0.41*log(yplus(i))+5.5;
    end
end
%semilogx(DelAlamo550Profiles.yplus,DelAlamo550Profiles.uplus,'+k','Markersize',2);
semilogx(yplus,theory,'k','linewidth',plotSpecs.linewidth);

hold
for plotNum = runsCup
    WCplottingFunctions(dataWC(plotNum).umean,y,'timeAveProfilePlot',plotNum,params,plotSpecs,0,0,1)
end

xlab = xlabel('\textbf{$z^{+}$}');
ylab = ylabel('\textbf{$u_{c}^{+}$}');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize)

leg = legend('log-law','run 200WC','run 350WC',...
    'run 500WC','location','northwest');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize)
if savePlot
    print('./Figures/timeAvesemilogU','-depsc')
    print('./jpegs/timeAvesemilogU','-djpeg','-r600')
end

%%%%%Sediment time series

fig2 = figure;
fig2.PaperUnits = 'centimeters';
fig2.PaperPosition = plotSpecs.fullWidthFigSize;
set(gca,'box','on')
hold

for plotNum = runsAll
    xnorm = dataWC(plotNum).Twave/(dataWC(plotNum).dt*dataWC(plotNum).nsavePro);
    ynorm = mean(dataWC(plotNum).CsedTotal1);
    %WCplottingFunctions(1:length(dataWC(plotNum).CsedTotal1),dataWC(plotNum).CsedTotal1*1000,'timePlot',plotNum,params,plotSpecs,temp,0,0)
    WCplottingFunctions(1:length(dataWC(plotNum).CsedTotal1),(dataWC(plotNum).CsedTotal1-ynorm)*100,'periodAvePlot',plotNum,params,plotSpecs,xnorm,ynorm,0)
end


ylab = ylabel('$\frac{\langle \overline{M_{T}} \rangle_{T}-\langle \overline{M_{T}} \rangle_{all}}{\langle \overline{M_{T}} \rangle_{all}}$ (\%)');
%set(gca,'ylim',[min(y/(params.molecular_viscosity/(Re_max*params.molecular_viscosity/H))) max(y/(params.molecular_viscosity/(Re_max*params.molecular_viscosity/H)))+100],'xlim',[-0.25 0.25])
xlab = xlabel('$t/T$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize)

% h(1) = plot(0,0,'color',plotSpecs.C(1,:));%, 'visible', 'off');
% h(2) = plot(0,0,'color',plotSpecs.C(3,:));%, 'visible', 'off');
% h(3) = plot(0,0,'color',plotSpecs.C(5,:));%, 'visible', 'off');
% leg = legend(h,'$Re_{\tau} = 200$ + waves','$Re_{\tau} = 350$ + waves',...
%     '$Re_{\tau} = 500$ + waves','location','northwest');
% set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize)

if savePlot
    print('./Figures/TotalmassChange','-depsc')
    print('./jpegs/TotalmassChange','-djpeg','-r600')
end

%%%%%Total Kinetic Energy time series

fig222 = figure;
fig222.PaperUnits = 'centimeters';
fig222.PaperPosition = plotSpecs.fullWidthFigSize;
set(gca,'box','on')
hold

for plotNum = runsAll
    xnorm = dataWC(plotNum).Twave/(dataWC(plotNum).dt*dataWC(plotNum).nsavePro);
    ynorm = mean(dataWC(plotNum).totalKinetic);
    WCplottingFunctions(1:length(dataWC(plotNum).totalKinetic),(dataWC(plotNum).totalKinetic-ynorm)*100,'periodAvePlot',plotNum,params,plotSpecs,xnorm,ynorm,0)
end


ylab = ylabel('$\frac{\langle \overline{K} \rangle_{T}-\langle \overline{K} \rangle_{all}}{\langle \overline{K} \rangle_{all}}$ (\%)');
xlab = xlabel('$t/T$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize)

h(1) = plot(0,0,'color',plotSpecs.C(1,:));%, 'visible', 'off');
h(2) = plot(0,0,'color',plotSpecs.C(3,:));%, 'visible', 'off');
h(3) = plot(0,0,'color',plotSpecs.C(5,:));%, 'visible', 'off');
leg = legend(h,'$Re_{\tau} = 200$ + waves','$Re_{\tau} = 350$ + waves',...
    '$Re_{\tau} = 500$ + waves','location','northwest');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize)

if savePlot
    print('./Figures/TotalKineticEnergyChange','-depsc')
    print('./jpegs/TotalKineticEnergyChange','-djpeg','-r600')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%U_depthAve%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig3 = figure;
fig3.PaperUnits = 'centimeters';
fig3.PaperPosition = plotSpecs.fullWidthFigSize;
set(gca,'box','on')
hold

for plotNum = runsAll
    ustar_temp = sqrt(dataWC(plotNum).dpdxSteady*H/params.rhoWater);
    ynorm = ustar_temp/kappa*(log(H/(params.molecular_viscosity/(9*ustar_temp)))+params.molecular_viscosity/(9*ustar_temp)/H-1);
    xnorm = dataWC(plotNum).Twave/(dataWC(plotNum).dt*dataWC(plotNum).nsavePro);
    WCplottingFunctions(1:length(dataWC(plotNum).uDepthAve),(dataWC(plotNum).uDepthAve -ynorm)*100,'periodAvePlot',plotNum,params,plotSpecs,xnorm,ynorm,0)
end


ylab = ylabel('$\frac{\langle \overline{u} \rangle_{T} - u_{0}}{u_{0}}$');
xlab = xlabel('$t/T$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize)

h(1) = plot(0,0,'color',plotSpecs.C(1,:));%, 'visible', 'off');
h(2) = plot(0,0,'color',plotSpecs.C(3,:));%, 'visible', 'off');
h(3) = plot(0,0,'color',plotSpecs.C(5,:));%, 'visible', 'off');
leg = legend(h,'$Re_{\tau} = 200$ + waves','$Re_{\tau} = 350$ + waves',...
    '$Re_{\tau} = 500$ + waves','location','southeast');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize)


if savePlot
    print('./Figures/depthAveVel','-depsc')
    print('./jpegs/depthAveVel','-djpeg','-r600')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Velocity Phase Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig3 = figure;
fig3.PaperUnits = 'centimeters';
fig3.PaperPosition = plotSpecs.quaterFigSize;
set(gca,'box','on')
hold

for phase = 1:5:60
    for plotNum = runsCup
        xnorm = dataWC(plotNum).dpdxWave/(params.rhoWater*2*pi/dataWC(plotNum).Twave);
        ynorm = H;
        WCplottingFunctions(dataWC(plotNum).umean,y,'phaseAvePlot',plotNum,params,plotSpecs,xnorm,ynorm,phase)
    end
end


ylab = ylabel('$\frac{z}{H}$');
%set(gca,'ylim',[min(y/(params.molecular_viscosity/(Re_max*params.molecular_viscosity/H))) max(y/(params.molecular_viscosity/(Re_max*params.molecular_viscosity/H)))+100],'xlim',[-0.25 0.25])
xlab = xlabel('$\frac{\widetilde{u}}{u_{b}}$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize)

% h(1) = plot(0,0,'color',plotSpecs.C(1,:));%, 'visible', 'off');
% h(2) = plot(0,0,'color',plotSpecs.C(3,:));%, 'visible', 'off');
% h(3) = plot(0,0,'color',plotSpecs.C(5,:));%, 'visible', 'off');
% leg = legend(h,'$Re_{\tau} = 200$ + waves','$Re_{\tau} = 350$ + waves',...
%     '$Re_{\tau} = 500$ + waves','location','northeast');
% set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize)

if savePlot
    print('./Figures/phaseAveVel','-depsc')
    print('./jpegs/phaseAveVel','-djpeg','-r600')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Wave component of phase averaged velcoity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig4 = figure;
fig4.PaperUnits = 'centimeters';
fig4.PaperPosition = plotSpecs.quaterFigSize;
set(gca,'box','on')
hold
xlim = [-1.2 1.2];
ylim = [0 0.15];
        count = 0;
for phase = 1:5:60
    for plotNum = runsCup    
        xnorm = dataWC(plotNum).dpdxWave/(params.rhoWater*2*pi/dataWC(plotNum).Twave);
        ynorm = H;
        WCplottingFunctions(dataWC(plotNum).umean,y,'WaveComPhase',plotNum,params,plotSpecs,xnorm,ynorm,phase)
    end
        params.delta = sqrt(2*params.molecular_viscosity/(2*pi/dataWC(plotNum).Twave));
        WCplottingFunctions(0,y,'WaveTheo',0,params,0,1,H,count)
        count = count+pi/6;
end

ylab = ylabel('$\frac{z}{H}$');
xlab = xlabel('$\frac{\widetilde{u}_{w}}{u_{b}}$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize,'ylim',ylim,'xlim',xlim)

if savePlot
    print('./Figures/phaseAveWaveVel','-depsc')
    print('./jpegs/phaseAveWaveVel','-djpeg','-r600')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Planform and wave component of phase averaged velcoity side by side%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig34 = figure;
fig34.PaperUnits = 'centimeters';
fig34.PaperPosition = plotSpecs.fullWidthHalfPage;
clf
ha = tight_subplot(1,2,[0.02 .05],[.15 .04],[.08 .04]);
count = 1;
ylim = [0 8];
xlim1 = [-2 2];
xlim2 = [-2 2];

axes(ha(1))
set(gca,'box','on')
hold on
for phase = 1:5:60
    for plotNum = 7
        xnorm = dataWC(plotNum).dpdxWave/(params.rhoWater*2*pi/dataWC(plotNum).Twave);
        ynorm = delta;
        WCplottingFunctions(dataWC(plotNum).umean,y,'phaseAvePlot',plotNum,params,plotSpecs,xnorm,ynorm,phase)
    end
end

for phase = 1:5:60
    for plotNum = [runsCup, 7]
        xnorm = dataWC(plotNum).dpdxWave/(params.rhoWater*2*pi/dataWC(plotNum).Twave);
        ynorm = delta;
        WCplottingFunctions(dataWC(plotNum).umean,y,'phaseAvePlot',plotNum,params,plotSpecs,xnorm,ynorm,phase)
    end
end


ylab = ylabel('$z/\Delta$');
xlab = xlabel('$\langle \widetilde{u} \rangle_{\hbox{p}}/u_{b}$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize,'ylim',ylim,'xlim',xlim1)


axes(ha(2))
set(gca,'box','on')
hold on
for phase = 1:5:60
    for plotNum = 7   
        xnorm = dataWC(plotNum).dpdxWave/(params.rhoWater*2*pi/dataWC(plotNum).Twave);
        ynorm = delta;
        WCplottingFunctions(dataWC(plotNum).umean,y,'WaveComPhase',plotNum,params,plotSpecs,xnorm,ynorm,phase)
    end
end


for phase = 1:5:60
    for plotNum = runsCup   
        xnorm = dataWC(plotNum).dpdxWave/(params.rhoWater*2*pi/dataWC(plotNum).Twave);
        ynorm = delta;
        WCplottingFunctions(dataWC(plotNum).umean,y,'WaveComPhase',plotNum,params,plotSpecs,xnorm,ynorm,phase)
    end
end

count = 0;
for phase = 1:5:60
        params.delta = sqrt(2*params.molecular_viscosity/(2*pi/dataWC(plotNum).Twave));
        WCplottingFunctions(0,y,'WaveTheo',0,params,0,1,delta,count)
        count = count+pi/6;
end

xlab2 = xlabel('$\langle u_{\hbox{w}} \rangle_{\hbox{p}}/u_{b}$');
set(xlab2,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize,'ylim',ylim,'xlim',xlim2,'yticklabel',[])

h(1) = plot(0,0,'xk','markersize', 3);%, 'visible', 'off');
h(2) = plot(0,0,'color',plotSpecs.C(7,:));%, 'visible', 'off');
h(3) = plot(0,0,'color',plotSpecs.C(1,:));%, 'visible', 'off');
h(4) = plot(0,0,'color',plotSpecs.C(3,:));%, 'visible', 'off');
h(5) = plot(0,0,'color',plotSpecs.C(5,:));%, 'visible', 'off');
leg = legend(h,'Stokes solution','0WC','200WC','350WC','500WC');

set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.76779741559709 0.799166667745248 0.155297797066825 0.114285712582724])
a1 = annotation(fig34,'textbox',...
 [0.0930224084000881 0.903538216688553 0.0497025115149362 0.0409523810659136],...  
    'String','(a)',...
    'Interpreter','latex',...
    'FontSize',11,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);
set(a1,'interpreter','Latex')

a2 = annotation(fig34,'textbox',...
 [0.559093836971517 0.905919169069504 0.0497025115149362 0.0409523810659136],...  
    'String','(b)',...
    'Interpreter','latex',...
    'FontSize',11,...
    'FitBoxToText','on',...
     'LineStyle','none',...
    'Margin',2,...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);
set(a2,'interpreter','Latex')


if savePlot
    print('./Figures/planFromAndWaveU','-depsc')
    print('./jpegs/planFromAndWaveU','-djpeg','-r600')
end


% %%%%%%%%%%%%%%%%%%%Phase Panel Plot Vel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig5 = figure;
fig5.PaperUnits = 'centimeters';
fig5.PaperPosition = plotSpecs.fullPageFigSize;
set(gca,'box','on')
clf
ha = tight_subplot(4,3,[0.03 .02],[.1 .04],[.12 .08]);
count = 1;
ylim = [0 1];
xlim = [-1 2.5];

for phase = 1:5:60
    axes(ha(count))
    
    if count == 1 || count == 4 || count == 7 || count == 10
        ylab = ylabel('$\frac{z}{H}$');
        
        set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
        set(gca,'ylim',ylim)
        
    else
        set(gca,'yticklabel',[],'ylim',ylim)
    end
    
    if count == 10 || count == 11 || count == 12
        if count == 11
            xlab = xlabel('$\frac{\widetilde{u}}{u_{b}}$');
            set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
        end
        set(gca,'xlim',xlim)
    else
        set(gca,'xticklabel',[],'xlim',xlim)
    end
    
    hold
    for plotNum = runsCup
        xnorm = dataWC(plotNum).dpdxWave/(params.rhoWater*2*pi/dataWC(plotNum).Twave);
        ynorm = H;
        WCplottingFunctions(dataWC(plotNum).umean,y,'phaseAvePlot',plotNum,params,plotSpecs,xnorm,ynorm,phase)
    end
    set(gca,'FontSize',plotSpecs.axesfontsize)
    
    p = get(gca, 'Position');
    if count==1
        ax = axes('Parent', gcf, 'Position', [0.108001498186406 0.912242686890574 0.101501741554415 0.0421835240872219]);
        p2 = plot(ax, 1:12,cos([1:12]*pi/6),'k');
        set(ax, 'Xlim', [0 13]);
        set(ax,'XTick', [],'YTick', []);
        set(gca,'Visible','off')
        hold
        plot(1:12,cos([1:12]*pi/6),'k',1:12,cos([1:12]*pi/6),'*k','MarkerSize',2)
    end
    count = count+1;
end
if savePlot
    print('./Figures/phaseAveVelPanel','-depsc')
    print('./jpegs/phaseAveVelPanel','-djpeg','-r600')
end