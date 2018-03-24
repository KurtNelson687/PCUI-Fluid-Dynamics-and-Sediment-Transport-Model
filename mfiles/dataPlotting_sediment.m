% Script: dataPlotting_sediment
%
% Author: Kurt Nelson
%
% Purpose: This script creates 5 plots related to the sediment concentration and
% The script requires processed PCUI data from "dataExtractor.m".
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
runsAll = 1:6; % Runs to plot for wave and current simulations for stratified and unstratified runs
aveStart = 1; % Index to start phase averaging at
num_steps_period = 12*5; % Number of times data is recorded per wave cycle
phases = (0:num_steps_period)/num_steps_period*2*pi; % Number of phases to plot over
getDataAndGridInfo % Load PCUI data, plotting variables, grid data, and Moser and del Alamo data
plotSpecs.lineType = {'-','--',':'}; % Line types specific to this script
delta = sqrt(2*params.molecular_viscosity/(2*pi/3)); % Stokes layer thickness
kappa = 0.41; % von Karman constant
buffTop = 30; % z+ value chosen for top of the buffer layer
% x labels
xticklabels = {'0','\pi/4', '\pi/2','3\pi/4', '\pi','5\pi/4', '3\pi/2','7\pi/4', '2 \pi'};
xtick = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4,2*pi];

%% Bottom stress and bed concentrations
fig33 = figure;
fig33.PaperUnits = 'centimeters';
fig33.PaperPosition =  plotSpecs.quaterFigSize;
set(gca,'box','on')
xlim = [0 2*pi];

for plotNum = [runsCup, 7]
    temp = params.rhoWater*params.molecular_viscosity*dataWC(plotNum).umean(1,:)./y(1);
    phaseAveragedCsed =postProcess(dataWC(plotNum).Csed,x,y,z,'phaseAverage');
    phaseAveragedShear =postProcess(temp,x,y,z,'phaseAverage');
    
    if plotNum == runsCup(1)
        [AX,H1,H2] =plotyy(phases,[phaseAveragedCsed(1,:), phaseAveragedCsed(1,1)]*1000,...
            phases,[phaseAveragedShear, phaseAveragedShear(1)]);
        hold(AX(1))
        hold(AX(2))
    else
        plot(AX(1),phases,[phaseAveragedCsed(1,:), phaseAveragedCsed(1,1)]*1000,...
            'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
        plot(AX(2),phases, [phaseAveragedShear, phaseAveragedShear(1)],...
            '--','color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth);
    end
end
set(H1,'color',plotSpecs.C(1,:),'linewidth',plotSpecs.linewidth);
set(H2,'color',plotSpecs.C(1,:),'linewidth',plotSpecs.linewidth,'LineStyle','--');
set(AX,{'ycolor'},{'k';'k'})
set(AX(1),'ylim',[0 600],'ytick',[0,150,300,450,600],'FontSize',plotSpecs.axesfontsize,'xlim',xlim)
set(AX(2),'FontSize',plotSpecs.axesfontsize,'xlim',xlim)
set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel',xticklabels)

plot(AX(2),xlim,[-dataWC(plotNum).tauCrit -dataWC(plotNum).tauCrit ],':k','linewidth',1);
plot(AX(2),xlim,[dataWC(plotNum).tauCrit  dataWC(plotNum).tauCrit],':k','linewidth',1);
plot(AX(2),xlim,[0  0],':k','linewidth',1);


set(gca,'FontSize',plotSpecs.axesfontsize)
ylabel(AX(1),'$\langle \widetilde{C}_{b} \rangle_{\mathrm{p}}$','interpreter','Latex','FontSize',plotSpecs.labfontsize) % left y-axis
ylabel(AX(2),'$\langle \tau_{\mathrm{b}} \rangle_{\mathrm{p}}$','interpreter','Latex','FontSize',plotSpecs.labfontsize) % left y-axis
xlab = xlabel('\textbf{$ \theta $}');
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)

h(1) = plot(0,0,'color',plotSpecs.C(7,:));%, 'visible', 'off');
h(2) = plot(0,0,'color',plotSpecs.C(1,:));%, 'visible', 'off');
h(3) = plot(0,0,'color',plotSpecs.C(3,:));%, 'visible', 'off');
h(4) = plot(0,0,'color',plotSpecs.C(5,:));%, 'visible', 'off');
leg = legend(h,'0W','200WC','350WC','500WC');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize)

if savePlot
    print('./Figures/bottomStressAndBedSed','-depsc')
    print('./jpegs/bottomStressAndBedSed','-djpeg','-r600')
end


%% Bed Stress and Entrainment
fig333 = figure;
fig333.PaperUnits = 'centimeters';
fig333.PaperPosition =  [0 0 14 10];
clf
ha = tight_subplot(3,1,[0.03 .05],[.1 .125],[.1 .1]); %axes for plots
xlim = [0 2*pi];

axes('Position',[0.1 0.897619047619048 0.8 0.0895238095238177]) % position of top plot
set(gca,'box','on')
plot(phases,sin(phases),'-',phases,cos(phases),'--')
set(gca,'xtick',xtick,'xticklabel',[],'xlim',xlim,'ylim',[-1.1 1.1])
ylab = ylabel('$\frac{u_\mathrm{\infty}}{u_\mathrm{b}}$ and $\frac{S_\mathrm{W}}{u_\mathrm{b} \omega}$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel','')

leg = legend('$u_\mathrm{\infty}$','$S_\mathrm{W}$');

set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.139282370465142 0.909047619217921 0.175006687641144 0.0307142853736877],...
    'Orientation','horizontal')

axes(ha(1))
set(gca,'box','on')
hold

for plotNum = [7, runsCup]
    temp =    params.rhoWater*params.molecular_viscosity*dataWC(plotNum).umean(1,:)./y(1);
    phaseAveragedShear =postProcess(temp,x,y,z,'phaseAverage');
    plot(phases,[phaseAveragedShear, phaseAveragedShear(1)]/dataWC(plotNum).tauCrit,...
        'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
end
plot(xlim,[-1 -1],':k','linewidth',1);
plot(xlim,[1 1],':k','linewidth',1);
ylab = ylabel('$\frac{\langle \tau |_{z=0} \rangle_{\mathrm{p}}} {\tau_{crit}}$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel','','xlim',xlim)
leg = legend('0W','200WC','350WC','500WC');

set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.394270631245197 0.831803522507362 0.416716572216579 0.0307142853736877],...
    'Orientation','horizontal')

axes(ha(2))
set(gca,'box','on')
hold
for plotNum = [runsCup, 7]
    temp =    params.rhoWater*params.molecular_viscosity*dataWC(plotNum).umean(1,:)./y(1);
    phaseAveragedCsed =postProcess(dataWC(plotNum).Csed,x,y,z,'phaseAverage');
    phaseAveragedShear =postProcess(temp,x,y,z,'phaseAverage');
    index = abs(phaseAveragedShear)>=dataWC(plotNum).tauCrit;
    dataWC(plotNum).E = zeros(size(phaseAveragedShear));
    dataWC(plotNum).E(index) = dataWC(plotNum).Ased*(abs(phaseAveragedShear(index))-dataWC(plotNum).tauCrit);
    D = phaseAveragedCsed(1,:)*dataWC(plotNum).ws;
    Etot = (dataWC(plotNum).E-D)/(dataWC(plotNum).Ased*dataWC(plotNum).tauCrit);
    plot(phases,[Etot, Etot(1)],...
        'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
end
ylab = ylabel('$\frac{\langle \widetilde{E}_{\mathrm{net}} \rangle_{\mathrm{p}}}{M \tau_{crit}}$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel','','xlim',xlim)

axes(ha(3))
set(gca,'box','on')
hold
for plotNum = [runsCup, 7]
    
    if plotNum == 7
        phaseAveragedCsed =postProcess(dataWC(plotNum).Csed,x,y,z,'phaseAverage');
    else
        phaseAveragedCsed =postProcess(dataWC(plotNum).Csed(:,6361:end),x,y,z,'phaseAverage');
    end
    phaseAveragedCtot =postProcess( phaseAveragedCsed,x,y,z,'depthAverage');
    
    plot(phases,[phaseAveragedCtot,phaseAveragedCtot(1)]/mean(phaseAveragedCtot),...
        'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
end
ylab = ylabel('$\frac{\langle \overline{C} \rangle_{\mathrm{p}}}{\langle \overline{C} \rangle}$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel',xticklabels,'xlim',xlim)

xlab = xlabel('\textbf{$ \theta $}');
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)

addText('(a)', fig333,[0.102831319877079 0.835238095408395 0.0328373602458409 0.0295238091832115]); % text for panel label
addText('(b)', fig333,[0.102136863606317 0.566190476360782 0.0342262727873666 0.0295238091832115]); % text for panel label
addText('(c)', fig333,[0.0992598913056509 0.306666666836977 0.0328373602458409 0.0295238091832116]); % text for panel label



if savePlot
    print('./Figures/bottomStressEntrainment','-depsc')
    print('./jpegs/bottomStressEntrainment','-djpeg','-r600')
end

%% Time average sediment concentration
fig4 = figure;
fig4.PaperUnits = 'centimeters';
fig4.PaperPosition = plotSpecs.quaterFigSize;
set(gca,'box','on')
hold
xlim = [0 100];
ylim = [0 105];
for plotNum = [runsCup, 7]
    WCplottingFunctions(dataWC(plotNum).Csed(15:end,:)*1000,y(15:end),'timeAveProfilePlot',plotNum,params,plotSpecs,1,delta,0)
end

% plot theoretical Rouse profile taking Ca as the period average
% sediment concentration at the height given by matchIndex
for i = [7,2,4,6]
    
    % compute period averaged sediment concentration
    [~,endAve] =size(dataWC(i).Csed);
    endAve = floor(endAve/num_steps_period)*num_steps_period;
    aveIndex = 1:endAve;
    
    if i == 7
        ytheo = y(15:5:end); % height vector for theoretical solution
        Ca = mean(dataWC(i).Csed(1,aveIndex),2)'*1000; % time averaged concentration at the bed
        Ctheo = Ca*exp(-dataWC(i).ws/dataWC(i).ak*ytheo); % wave sediment profile
        plot(Ctheo,ytheo/delta,'ok','markersize',3)
    else % compute theoretical Rouse profile for wave and current runs
        u_star_plot = sqrt(dataWC(i).dpdxSteady*H/params.rhoWater); % u_star_plot
        matchIndex(i)=find(y/(params.molecular_viscosity/u_star_plot)>buffTop,1,'first'); % index of top of buffer layer
        ytheo = y(matchIndex(i)):0.005:H; % height vector for theoretical solution
        Ro = dataWC(i).ws/(kappa*u_star_plot); % Rouse number
        Ca = mean(dataWC(i).Csed(:,aveIndex),2)'; % Sediment concentration at the top of the buffer layer
        Ca = Ca(matchIndex(i))*1000;
        Ctheo = Ca*((H-ytheo)./ytheo*ytheo(1)./(H-ytheo(1))).^Ro; % Theoretical Rouse Profile
        plot(Ctheo,ytheo/delta,'xk','markersize',3)
        CsedBedToBuffer = mean(dataWC(i).Csed(1,aveIndex),2)'*1000/Ca
    end
end

ylab = ylabel('$z/\Delta$');

xlab = xlabel('$C_{\mathrm{c}}$ ($\frac{mg}{L}$)');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize,'xlim',xlim,'ylim',ylim)

leg = legend('0W','200WC','350WC','500WC','Theoretical wave solution','Rouse profile');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize)

if savePlot
    print('./Figures/timeAveSed','-depsc')
    print('./jpegs/timeAveSed','-djpeg','-r600')
end

%% Difference in time average sediment concentration
fig44 = figure;
fig44.PaperUnits = 'centimeters';
fig44.PaperPosition = plotSpecs.quaterFigSize;
set(gca,'box','on')
hold
num_steps_period = 12*5;

for plotNum = 1:floor(length(runsAll)/2)
    [~,endAve] =size(dataWC(plotNum*2).Csed);
    endAve = floor(endAve/num_steps_period)*num_steps_period;
    aveIndex = 60*70+1:endAve;
    
    xplot = mean(dataWC(plotNum*2-1).Csed(:,aveIndex),2)-mean(dataWC(plotNum*2).Csed(:,aveIndex),2);
    depthAvePercentDiff(plotNum) = y*xplot/(y*mean(dataWC(plotNum*2-1).Csed(:,aveIndex),2))
    xplot =xplot*100;
    xplotNorm = mean(dataWC(plotNum*2).Csed(:,aveIndex),2);
    plot(xplot./xplotNorm,y/delta,...
        'color',plotSpecs.C(plotNum*2,:),'linewidth',0.5)
end

ylab = ylabel('$z/\Delta$');
xlab = xlabel('$\frac{\langle \widetilde{C}_{no \, \rho} \rangle_{T} -\langle \widetilde{C}_{\rho} \rangle_{T}}{\langle \widetilde{C}_{\rho} \rangle_{T}}$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize)

if savePlot
    print('./Figures/diffTimeAveSed','-depsc')
    print('./jpegs/diffTimeAveSed','-djpeg','-r600')
end


%%%%%%%Time average sediment concentration and difference all%%%%%%%
fig444 = figure;
fig444.PaperUnits = 'centimeters';
fig444.PaperPosition =  plotSpecs.fullWidthFigSize;
set(gca,'box','on')
clf
ha = tight_subplot(1,2,[0.02 .05],[.17 .04],[.08 .04]);
count = 1;
ylim1 = [0 1];
xlim1 = [0 600];
axes(ha(1))
hold on
for plotNum = runsAll
    WCplottingFunctions(dataWC(plotNum).Csed*1000,y,'timeAveProfilePlot',plotNum,params,plotSpecs,1,H,0)
end

ylab = ylabel('$\frac{z}{H}$');
xlab = xlabel('$\langle \widetilde{C} \rangle_{T}$ ($\frac{mg}{L}$)');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize)


axes(ha(2))
hold on
for plotNum = 1:floor(length(runsAll)/2)
    [~,endAve] =size(dataWC(plotNum*2).Csed);
    endAve = floor(endAve/num_steps_period)*num_steps_period;
    aveIndex = 1:endAve;
    
    xplot = mean(dataWC(plotNum*2-1).Csed(:,aveIndex),2)-mean(dataWC(plotNum*2).Csed(:,aveIndex),2);
    xplot =xplot;
    xplotNorm = mean(dataWC(plotNum*2).Csed(:,aveIndex),2);
    plot(xplot./xplotNorm,y/H,...
        'color',plotSpecs.C(plotNum*2,:),'linewidth',plotSpecs.linewidth)
end

xlab = xlabel('$\frac{\langle \widetilde{C}_{no \, \rho} \rangle_{T} -\langle \widetilde{C}_{\rho} \rangle_{T}}{\langle \widetilde{C}_{\rho} \rangle_{T}}$');
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize,'yticklabel',[],'xlim',[0 0.80])

if savePlot
    print('./Figures/timeAveSedandSedDiff','-depsc')
    print('./jpegs/timeAveSedandSedDiff','-djpeg','-r600')
end


%%%%%%%Time average vertical turbulent sediment flux%%%%%%%
fig48 = figure;
fig48.PaperUnits = 'centimeters';
fig48.PaperPosition = plotSpecs.quaterFigSize;
set(gca,'box','on')
hold

for plotNum = runsAll
    %xnorm = sqrt(dataWC(plotNum).dpdxSteady*H/params.rhoWater)*dataWC(plotNum).rhoSed;
    xnorm = 1;
    ynorm = delta;
    WCplottingFunctions(-1*dataWC(plotNum).vCsed,y,'timeAveProfilePlot',plotNum,params,plotSpecs,xnorm,ynorm,0)
end

ylab = ylabel('$z/\Delta$');
xlab = xlabel('$(\widetilde{w''C''})_{\mathrm{c}}$ ($\frac{kg}{m^{2} s})$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
ylim = [0 105];
xlim = [0 3.5*10^-5];
set(gca,'FontSize',plotSpecs.axesfontsize,'Xlim',xlim,'Ylim',ylim);

h(1) = plot(0,0,'color',plotSpecs.C(1,:));%, 'visible', 'off');
h(2) = plot(0,0,'--','color',plotSpecs.C(1,:));%, 'visible', 'off');
h(3) = plot(0,0,'color',plotSpecs.C(3,:));%, 'visible', 'off');
h(4) = plot(0,0,'--','color',plotSpecs.C(3,:));%, 'visible', 'off');
h(5) = plot(0,0,'color',plotSpecs.C(5,:));%, 'visible', 'off');
h(6) = plot(0,0,'--','color',plotSpecs.C(5,:));%, 'visible', 'off');

leg = legend(h,'200WCa','200WCb','350WCa','350WCb','500WCa','500WCb');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize)

if savePlot
    print('./Figures/timeAvewCsed','-depsc')
    print('./jpegs/timeAvewCsed','-djpeg','-r600')
end

%%%%%%%Difference in time average vertical turbulent sediment flux%%%%%%%
fig49 = figure;
fig49.PaperUnits = 'centimeters';
fig49.PaperPosition = plotSpecs.quaterFigSize;
set(gca,'box','on')
hold
num_steps_period = 12*5;

for plotNum = 1:floor(length(runsAll)/2)
    [~,endAve] =size(dataWC(plotNum*2).Csed);
    endAve = floor(endAve/num_steps_period)*num_steps_period;
    aveIndex = 1:endAve;
    
    xplot = -1*mean(dataWC(plotNum*2-1).vCsed(:,aveIndex),2)+mean(dataWC(plotNum*2).vCsed(:,aveIndex),2);
    xplot =xplot*100;
    %xplotNorm = -1*mean(dataWC(plotNum*2).vCsed(:,aveIndex),2);
    xplotNorm = 1;
    plot(xplot./xplotNorm,y/H,...
        'color',plotSpecs.C(plotNum*2,:),'linewidth',0.5)
end

ylab = ylabel('$\frac{z}{H}$');
xlab = xlabel('$(\langle \widetilde{w''C''}_{no \, \rho} \rangle_{T}- \langle \widetilde{w''C''}_{\rho} \rangle_{T})$');

%xlab = xlabel('$\frac{\widetilde{w''C''}_{no \, \rho}-\widetilde{w''C''}_{\rho}}{\widetilde{w''C''}_{\rho}}$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
ylim = [0 1];
xlim = [0 1.5*10^-4];
%xlim = [0 100];
set(gca,'FontSize',plotSpecs.axesfontsize,'Xlim',xlim,'Ylim',ylim);

% h(1) = plot(0,0,'color',plotSpecs.C(1,:));%, 'visible', 'off');
% h(2) = plot(0,0,'color',plotSpecs.C(3,:));%, 'visible', 'off');
% h(3) = plot(0,0,'color',plotSpecs.C(5,:));%, 'visible', 'off');
% leg = legend(h,'$Re_{\tau} = 200$ + waves','$Re_{\tau} = 350$ + waves',...
%     '$Re_{\tau} = 500$ + waves','location','northeast');
% set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize)

if savePlot
    print('./Figures/diffTimeAveSedFlux','-depsc')
    print('./jpegs/diffTimeAveSedFlux','-djpeg','-r600')
end

% %%%%%%%%%%%%%%%%%%%Phase Panel Plot Sediment Concentration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig9 = figure;
fig9.PaperUnits = 'centimeters';
fig9.PaperPosition = plotSpecs.fullPageFigSize;
set(gca,'box','on')
clf
ha = tight_subplot(4,3,[0.03 .025],[.1 .04],[.12 .08]);
count = 1;
ylim = [0 5];
xlim = [0 600];

% ylim = [0 30];
% xlim = [-6 6];
for phase = 1:5:60
    axes(ha(count))
    set(gca,'box','on')
    
    if count == 1 || count == 4 || count == 7 || count == 10
        ylab = ylabel('$z/\Delta$');
        
        set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
        set(gca,'ylim',ylim)
        
    else
        set(gca,'yticklabel',[],'ylim',ylim)
    end
    
    if count == 10 || count == 11 || count == 12
        if count == 11
            xlab = xlabel('$\langle \widetilde{C}  \rangle_{\mathrm{p}}$ ($\frac{mg}{L}$)');
            %                xlab = xlabel('$(\langle \widetilde{C}  \rangle_{\mathrm{p}} - \langle \widetilde{C}  \rangle)/\langle \widetilde{C}  \rangle$');
            
            set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
        end
        %set(gca,'xlim',xlim)
    else
        %set(gca,'xticklabel',[],'xlim',xlim)
    end
    
    hold
    for plotNum = [runsCup, 7]
        temp = 1;
        %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %           [~,endAve] =size(dataWC(plotNum).Csed);
        %           endAve = floor(endAve/num_steps_period)*num_steps_period;
        %           aveIndex = 1:endAve;
        %           xplotNorm = mean(dataWC(plotNum).Csed(:,aveIndex),2);
        %
        %           WCplottingFunctions((dataWC(plotNum).Csed-xplotNorm)./xplotNorm*100,y,'phaseAvePlot',plotNum,params,plotSpecs,temp,delta,phase)
        % %%%%%%%%%%%%%%%%%%%%%%%To plot difference in sediment
        % %%%%%%%%%%%%%%%%%%%%%%%concentration%%%%%%%%%%%%%%%
        
        
        WCplottingFunctions(dataWC(plotNum).Csed*1000,y,'phaseAvePlot',plotNum,params,plotSpecs,temp,delta,phase)
    end
    set(gca,'FontSize',plotSpecs.axesfontsize)
    
    count = count+1;
end

a1 = annotation(fig9,'textbox',...
    [0.308421502847779 0.928197022251018 0.0466913870402745 0.0295238091832115],...
    'String','$\theta = 0$',...
    'Interpreter','latex',...
    'FontSize',7,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);

a2 = annotation(fig9,'textbox',...
    [0.549206752489743 0.922482736309666 0.100835173470633 0.0409523810659136],...
    'String','$\theta = \pi/6$',...
    'Interpreter','latex',...
    'FontSize',7,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);
a3 = annotation(fig9,'textbox',...
    [0.841241719094792 0.928197022251019 0.0667652402605329 0.0295238091832115],...
    'String','$\theta =  \pi/3$',...
    'Interpreter','latex',...
    'FontSize',7,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);

a4 = annotation(fig9,'textbox',...
    [0.291241719094793 0.706768450822446 0.0667652402605329 0.0295238091832115],...
    'String','$\theta = \pi/2$',...
    'Interpreter','latex',...
    'FontSize',7,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);

a5 = annotation(fig9,'textbox',...
    [0.562670290523365 0.706768450822447 0.0667652402605329 0.0295238091832115],...
    'String','$\theta = 2 \pi/3$',...
    'Interpreter','latex',...
    'FontSize',7,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);

a6 = annotation(fig9,'textbox',...
    [0.837670290523363 0.706768450822448 0.0667652402605329 0.0295238091832115],...
    'String','$\theta = 5 \pi/6$',...
    'Interpreter','latex',...
    'FontSize',7,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);

a7 = annotation(fig9,'textbox',...
    [0.309545319474192 0.485339879393875 0.048015182358878 0.0295238091832115],...
    'String','$\theta = \pi$',...
    'Interpreter','latex',...
    'FontSize',7,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);

a8 = annotation(fig9,'textbox',...
    [0.562670290523363 0.482958927012923 0.0667652402605329 0.0295238091832115],...
    'String','$\theta = 7 \pi/6$',...
    'Interpreter','latex',...
    'FontSize',7,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);

a9 = annotation(fig9,'textbox',...
    [0.837670290523363 0.485339879393876 0.0667652402605329 0.0295238091832115],...
    'String','$\theta = 4 \pi/3$',...
    'Interpreter','latex',...
    'FontSize',7,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);

a10 = annotation(fig9,'textbox',...
    [0.287670290523364 0.263911307965305 0.0667652402605329 0.0295238091832116],...
    'String','$\theta = 3 \pi/2$',...
    'Interpreter','latex',...
    'FontSize',7,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);

a11 = annotation(fig9,'textbox',...
    [0.562670290523348 0.261530355584353 0.0667652402605329 0.0295238091832116],...
    'String','$\theta = 5 \pi/3$',...
    'Interpreter','latex',...
    'FontSize',7,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);

a12 = annotation(fig9,'textbox',...
    [0.829188139151271 0.261530355584352 0.0730152572904315 0.0295238091832116],...
    'String','$\theta = 11 \pi/6$',...
    'Interpreter','latex',...
    'FontSize',7,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);

h(1) = plot(0,0,'color','k');%, 'visible', 'off');
h(2) = plot(0,0,'color',plotSpecs.C(2,:));%, 'visible', 'off');
h(3) = plot(0,0,'color',plotSpecs.C(4,:));%, 'visible', 'off');
h(4) = plot(0,0,'color',plotSpecs.C(6,:));%, 'visible', 'off');

leg = legend(h,'0W', '200WC','350WC','500WC');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.77321203576852 0.822975324885751 0.113492148263114 0.105879838126046])

if savePlot
    print('./Figures/phaseAveCsedPanel','-depsc')
    print('./jpegs/phaseAveCsedPanel','-djpeg','-r600')
end


%%%%%%%%%%%%%%%%%%%Phase Panel Plot vertical turbulent flux%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig10 = figure;
fig10.PaperUnits = 'centimeters';
fig10.PaperPosition = plotSpecs.fullPageFigSize;
set(gca,'box','on')
clf
ha = tight_subplot(4,3,[0.03 .02],[.1 .04],[.12 .08]);
count = 1;
ylim = [0 1];
xlim = [0 4*10^-5];
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
            xlab = xlabel('$-\widetilde{w''C''}$ ($\frac{kg}{m^{2} s}$)');
            set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
        end
        set(gca,'xlim',xlim)
    else
        set(gca,'xticklabel',[],'xlim',xlim)
    end
    
    hold
    for plotNum = runsCup
        temp = 1;
        WCplottingFunctions(-1*dataWC(plotNum).vCsed,y,'phaseAvePlot',plotNum,params,plotSpecs,temp,H,phase)
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
    print('./Figures/phaseAvevCsedPanel','-depsc')
    print('./jpegs/phaseAvevCsedPanel','-djpeg','-r600')
end
%% Functions
function addText(string, fig,position)
% adds text to figure specified by "string" to figure "fig" at position
% "position". Latex format is used.
annotation(fig,'textbox',...
    position,...
    'String',string,...
    'Interpreter','latex',...
    'FontSize',7,...
    'FitBoxToText','on',...
    'Margin',2,...
    'LineStyle','none',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center');
end



