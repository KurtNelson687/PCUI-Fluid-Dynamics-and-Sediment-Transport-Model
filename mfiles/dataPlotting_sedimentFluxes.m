% Script: dataPlotting_sedimentFluxes
%
% Author: Kurt Nelson
%
% Purpose: This script creates six plots related to sediment fluxes and
% divergence of sediment fluxes. The plots are 1) period-averged divergence
% of sediment fluxes, 2) period-averged sediment fluxes, 3) profiles of divergence of
% sediment fluxes form z = 0 to z = 5 delta, 4) profiles of divergence of
% sediment fluxes form z = 5 to 4 = delta, 5) fluxes integrated over layers
% with thickness delta from z = 0 to z = 5, and 6) bottom stress and fluxes
% integrated from the bed to the buffer layer, and from the buffer layer to
% the free-surface. The script requires processed PCUI data from "dataExtractor.m".
%%
close all; clear all; clear all
addpath('./Functions');

%% Flags
savePlot = true; %save plots or not

%% Variables
data_folder = '/Users/kurtnelson/Desktop/Writing/WaveAndCurrent_sedFocus/mfiles/Data'; %PCUI data folder
data_fileWC = 'dataWCfineSpinup'; %PCUI processed data file for currents only
data_fileC = 'dataCurrentOnly.mat'; %PCUI processed data file for currents only
runsCup = [2,4,6]; % Runs to plot
aveStart = 1; % Index to start phase averaging at
num_steps_period = 12*5; % Number of times data is recorded per wave cycle
phases = 1:1:60; % Index for phases
getDataAndGridInfo % Load PCUI data, plotting variables, grid data, and Moser and del Alamo data
plotSpecs.Cbudget = linspecer(3); % Line colors specific to this script
plotSpecs.lineType = {'-','--',':'}; % Line types specific to this script
delta = sqrt(2*params.molecular_viscosity/(2*pi/3)); % Stokes layer thickness
% x labels
xticklabels = {'0','\pi/4', '\pi/2','3\pi/4', '\pi','5\pi/4', '3\pi/2','7\pi/4', '2 \pi'};
xtick = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4,2*pi];
%% Assign coefficients for derivaties
%c1 is for first derivative
%c2 is for second derivative
[m] = length(y);
for j = 1:m
    if j == 1 %5-point: evaluation location, and 4 forward at bottom
        c1(j,:) = fdcoeffF(1,y(j),y(j:j+4));
        c2(j,:) = fdcoeffF(2,y(j),y(j:j+4));
    elseif j == 2 %5-point: 1 back, evaluation location, and 3 forward
        c1(j,:) = fdcoeffF(1,y(j),y(j-1:j+3));
        c2(j,:) = fdcoeffF(2,y(j),y(j-1:j+3));
    elseif j == m %5-point: evaluation location, and 4 back at top
        c1(j,:) = fdcoeffF(1,y(j),y(j-4:j));
        c2(j,:) = fdcoeffF(2,y(j),y(j-4:j));
    elseif j == m-1 %5-point: 3 back, evaluation location, and 1 forward
        c1(j,:) = fdcoeffF(1,y(j),y(j-3:j+1));
        c2(j,:) = fdcoeffF(2,y(j),y(j-3:j+1));
    else %5-point: 2 back, evaluation location, and 2 forward all internal points
        c1(j,:) = fdcoeffF(1,y(j),y(j-2:j+2));
        c2(j,:) = fdcoeffF(2,y(j),y(j-2:j+2));
    end
end
%% Compute current and wave sediment fluxes
dataWC = computeCurrentSedFlux(dataWC,c1,c2,runsCup,aveStart,num_steps_period);
dataWC = computeWaveSedFlux(dataWC,c1,c2,runsCup,aveStart,num_steps_period,phases);

%% Plot Divergence of current sediment fluxes
fig1 = figure;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = plotSpecs.quaterFigSize;
set(gca,'box','on')
hold

count = 1; %counter for line color
% Plot divergence of fluxes
for runNum = [2,4,6]
    plot(dataWC(runNum).steadyDiffusion(2:end),...
        y(2:end)/delta,plotSpecs.lineType{1},'color',plotSpecs.Cbudget(count,:),'linewidth',plotSpecs.linewidth)
    plot(dataWC(runNum).steadyDeposition(2:end),...
        y(2:end)/delta,plotSpecs.lineType{2},'color',plotSpecs.Cbudget(count,:),'linewidth',plotSpecs.linewidth)
    plot(dataWC(runNum).steadySedTurb(2:end),...
        y(2:end)/delta,plotSpecs.lineType{3},'color',plotSpecs.Cbudget(count,:),'linewidth',plotSpecs.linewidth)
    count=count+1;
end

myXlabel('Sediment Budget Terms (kg/(m$^3$ s))',plotSpecs.labfontsize)
myYlabel('$z/\Delta$',plotSpecs.labfontsize)

% legend
leg = legend('(200WC) $K \frac{\partial^2 C_{\hbox{c}}}{\partial z^2}$',...
    '(200WC) $\frac{\partial}{\partial z} (w_{s} C_{\hbox{c}})$',...
    '(200WC) -$\frac{\partial}{\partial z} (w''C''_{\hbox{c}})$',...
    '(350WC) $K \frac{\partial^2 C_{\hbox{c}}}{\partial z^2}$',...
    '(350WC) $\frac{\partial}{\partial z} (w_{s} C_{\hbox{c}})$',...
    '(350WC) -$\frac{\partial}{\partial z} (w''C''_{\hbox{c}})$',...
    '(500WC) $K \frac{\partial^2 C_{\hbox{c}}}{\partial z^2}$',...
    '(500WC) $\frac{\partial}{\partial z} (w_{s} C_{\hbox{c}})$',...
    '(500WC) -$\frac{\partial}{\partial z} (w''C''_{\hbox{c}})$');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.666003463949471 0.542133288440258 0.18753662109375 0.214816658269791])

%Plot over the bottom 5 delta
ylim = [0 5]; % limits for y axis
xlim = [-0.12 0.12]; % limits for x axis
set(gca,'FontSize',plotSpecs.axesfontsize,'xlim',xlim,'ylim',ylim)
if savePlot
    print('./Figures/periodSedBudget_bottom','-depsc')
    print('./jpegs/periodSedBudget_bottom','-djpeg','-r600')
end

%Plot over the bottom 100 delta
ylim = [1 100]; % y limits
xlim = [-4.5*10^-3 4.5*10^-3]; % x limits
set(gca,'FontSize',plotSpecs.axesfontsize,'xlim',xlim,'ylim',ylim)
if savePlot
    print('./Figures/periodSedBudget_top','-depsc')
    print('./jpegs/periodSedBudget_top','-djpeg','-r600')
end

%% Plot current sediment fluxes
fig2 = figure;
fig2.PaperUnits = 'centimeters';
fig2.PaperPosition = plotSpecs.quaterFigSize;
set(gca,'box','on')
hold

count = 1; % counter for line color
ylim = [0 100]; % y limits
xlim = [-4.5*10^-3 4.5*10^-3]; % x limits

% Plot fluxes
for runNum = [2,4,6]
    plot(-dataWC(runNum).Fdc(2:end),...
        y(2:end)/delta,plotSpecs.lineType{1},'color',plotSpecs.Cbudget(count,:),'linewidth',plotSpecs.linewidth)
    plot(-dataWC(runNum).Fwc(2:end),...
        y(2:end)/delta,plotSpecs.lineType{2},'color',plotSpecs.Cbudget(count,:),'linewidth',plotSpecs.linewidth)
    plot(-dataWC(runNum).Ftc(2:end),...
        y(2:end)/delta,plotSpecs.lineType{3},'color',plotSpecs.Cbudget(count,:),'linewidth',plotSpecs.linewidth)
    count=count+1;
end

myXlabel('Period-averaged sediment fluxes [kg/(m$^2$ s)]',plotSpecs.labfontsize)
myYlabel('$z/\Delta$',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize,'ylim',ylim) %'xlim',xlim,

%legend
leg = legend('$F_\mathrm{D,c}$ (200WC)',...
    '$F_\mathrm{s,c}$ (200WC)',...
    '$F_\mathrm{T,c}$ (200WC)',...
    '$F_\mathrm{D,c}$ (350WC)',...
    '$F_\mathrm{s,c}$ (350WC)',...
    '$F_\mathrm{T,c}$ (350WC)',...
    '$F_\mathrm{D,c}$ (500WC)',...
    '$F_\mathrm{s,c}$ (500WC)',...
    '$F_\mathrm{T,c}$ (500WC)');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.658462192331036 0.609541619107845 0.14904773575919 0.199047615982237])

if savePlot
    print('./Figures/periodSedFlux_top','-depsc')
    print('./jpegs/periodSedFlux_top','-djpeg','-r600')
end

%% Plot phase averaged divergence of sediment fluxes from 0 to 5 delta
fig3 = figure;
fig3.PaperUnits = 'centimeters';
fig3.PaperPosition = plotSpecs.fullPageFigSize;
set(gca,'box','on')
clf
ha = tight_subplot(4,3,[0.03 .025],[.1 .04],[.12 .08]);

count = 1; % counter for panel
ylim = [0 5]; % y limits
%xlim = [-4*10^-3 4*10^-3];
xlim = [-0.5 1]; % x limits
runNum =6; % run to plot

%plot phase panels for divergence of sediment fluxes
for phase = 1:5:60
    axes(ha(count))
    
    %y axis specification
    if count == 1 || count == 4 || count == 7 || count == 10
        myYlabel('$z/\Delta$',plotSpecs.labfontsize)
        set(gca,'ylim',ylim)
    else
        set(gca,'yticklabel',[],'ylim',ylim)
    end
    
    %x axis specification
    if count == 10 || count == 11 || count == 12
        if count == 11
            myXlabel('Sediment Budget Terms (kg/(m$^3$ s))',plotSpecs.labfontsize)
        end
        set(gca,'xlim',xlim)
    else
        set(gca,'xticklabel',[],'xlim',xlim)
    end
    
    hold
    
    plot(dataWC(runNum).steadySedTurb(1:end)+...
        dataWC(runNum).steadyDiffusion(1:end)+...
        dataWC(runNum).steadyDeposition(1:end)+...
        dataWC(runNum).waveSedTurb(1:end,phase)+...
        dataWC(runNum).waveDiffusion(1:end,phase)+...
        dataWC(runNum).waveDeposition(1:end,phase),...
        y(1:end)/delta,'k','linewidth',plotSpecs.linewidth)
    
    
    plot(dataWC(runNum).steadyDiffusion(1:end),...
        y(1:end)/delta,'--','color',plotSpecs.Cbudget(1,:),'linewidth',plotSpecs.linewidth)
    plot(dataWC(runNum).steadyDeposition(1:end),...
        y(1:end)/delta,'--','color',plotSpecs.Cbudget(2,:),'linewidth',plotSpecs.linewidth)
    plot(dataWC(runNum).steadySedTurb(1:end),...
        y(1:end)/delta,'--','color',plotSpecs.Cbudget(3,:),'linewidth',plotSpecs.linewidth)
    
    plot(dataWC(runNum).waveDiffusion(1:end,phase),...
        y(1:end)/delta,'color',plotSpecs.Cbudget(1,:),'linewidth',plotSpecs.linewidth)
    plot(dataWC(runNum).waveDeposition(1:end,phase),...
        y(1:end)/delta,'color',plotSpecs.Cbudget(2,:),'linewidth',plotSpecs.linewidth)
    plot(dataWC(runNum).waveSedTurb(1:end,phase),...
        y(1:end)/delta,'color',plotSpecs.Cbudget(3,:),'linewidth',plotSpecs.linewidth)
    
    set(gca,'FontSize',plotSpecs.axesfontsize)
    count = count+1;
end

% legend
clear h
h(1) = plot(0,0,'color','k');%, 'visible', 'off');
h(2) = plot(0,0,'--','color',plotSpecs.Cbudget(1,:));%, 'visible', 'off');
h(3) = plot(0,0,'--','color',plotSpecs.Cbudget(2,:));%, 'visible', 'off');
h(4) = plot(0,0,'--','color',plotSpecs.Cbudget(3,:));%, 'visible', 'off');
leg = legend(h,'$\frac{\partial \langle \widetilde{C} \rangle_{\hbox{p}}}{\partial t}$',...
    '$K \frac{\partial^2 C_{\hbox{c}}}{\partial z^2}$',...
    '$\frac{\partial}{\partial z} (w_{s} C_{\hbox{c}})$',...
    '-$\frac{\partial}{\partial z} (w''C'')_{\hbox{c}}$');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.836922805649892 0.860887412797838 0.120697937692915 0.105879838126046])

% Add text panel label for each panel
addText('$\theta = 0$', fig3,[0.308421502847779 0.928197022251018 0.0466913870402745 0.0295238091832115])
addText('$\theta = 1/6 \pi$', fig3,[0.549206752489743 0.922482736309666 0.100835173470633 0.0409523810659136])
addText('$\theta = 1/3 \pi$', fig3,[0.689456004809072 0.928197022251019 0.0667652402605329 0.0295238091832115])
addText('$\theta = 1/2 \pi$', fig3,[0.291241719094793 0.706768450822446 0.0667652402605329 0.0295238091832115])
addText('$\theta = 2/3 \pi$', fig3,[0.566241719094793 0.706768450822447 0.0667652402605329 0.0295238091832115])
addText('$\theta = 5/6 \pi$', fig3,[0.841241719094792 0.706768450822448 0.0667652402605329 0.0295238091832115])
addText('$\theta = \pi$', fig3,[0.309545319474192 0.485339879393875 0.048015182358878 0.0295238091832115])
addText('$\theta = 7/6 \pi$', fig3,[0.566241719094792 0.482958927012923 0.0667652402605329 0.0295238091832115])
addText('$\theta = 4/3 \pi$', fig3,[0.841241719094792 0.485339879393876 0.0667652402605329 0.0295238091832115])
addText('$\theta = 3/2 \pi$', fig3,[0.291241719094793 0.263911307965305 0.0667652402605329 0.0295238091832116])
addText('$\theta = 5/3 \pi$', fig3,[0.566241719094776 0.261530355584353 0.0667652402605329 0.0295238091832116])
addText('$\theta = 11/6 \pi$', fig3,[0.836330996294129 0.261530355584352 0.0730152572904315 0.0295238091832116])

% save profiles over the bottom 5 delta
if savePlot
    if runNum == 2
        print('./Figures/component_phaseSedBudget200WC_bottom','-depsc')
        print('./jpegs/component_phaseSedBudget200WC_bottom','-djpeg','-r600')
    elseif runNum == 4
        print('./Figures/component_phaseSedBudget350WC_bottom','-depsc')
        print('./jpegs/component_phaseSedBudget350WC_bottom','-djpeg','-r600')
    elseif runNum == 6
        print('./Figures/component_phaseSedBudget500WC_bottom','-depsc')
        print('./jpegs/component_phaseSedBudget500WC_bottom','-djpeg','-r600')
    end
end


%% Plot phase averaged divergence of sediment fluxes from 5 to 40 delta
fig4 = figure;
fig4.PaperUnits = 'centimeters';
fig4.PaperPosition = plotSpecs.fullPageFigSize;
set(gca,'box','on')
clf

ha = tight_subplot(4,3,[0.03 .025],[.1 .04],[.12 .08]);
xnorm = 10^-3; % normalization for the fluxes

count = 1; % counter for panels
ylim = [5 40]; % y limit
xlim = [-4*10^-3/xnorm 4*10^-3/xnorm]; % x limit

%plot phase panels for divergence of sediment fluxes
for phase = 1:5:60
    axes(ha(count))
    
    % x-label specification
    if count == 1 || count == 4 || count == 7 || count == 10
        myYlabel('$z/\Delta$',plotSpecs.labfontsize)
        set(gca,'ylim',ylim)
    else
        set(gca,'yticklabel',[],'ylim',ylim)
    end
    
    % y-label specification
    if count == 10 || count == 11 || count == 12
        if count == 11
            myXlabel('Sediment Budget Terms (kg/(m$^3$ s))',plotSpecs.labfontsize)
        end
        set(gca,'xlim',xlim)
    else
        set(gca,'xticklabel',[],'xlim',xlim)
    end
    
    hold
    
    plot((dataWC(runNum).steadySedTurb(1:end)+...
        dataWC(runNum).steadyDiffusion(1:end)+...
        dataWC(runNum).steadyDeposition(1:end)+...
        dataWC(runNum).waveSedTurb(1:end,phase)+...
        dataWC(runNum).waveDiffusion(1:end,phase)+...
        dataWC(runNum).waveDeposition(1:end,phase))/xnorm,...
        y(1:end)/delta,'k','linewidth',plotSpecs.linewidth)
    
    plot(dataWC(runNum).steadyDiffusion(1:end)/xnorm,...
        y(1:end)/delta,'--','color',plotSpecs.Cbudget(1,:),'linewidth',plotSpecs.linewidth)
    plot(dataWC(runNum).steadyDeposition(1:end)/xnorm,...
        y(1:end)/delta,'--','color',plotSpecs.Cbudget(2,:),'linewidth',plotSpecs.linewidth)
    plot(dataWC(runNum).steadySedTurb(1:end)/xnorm,...
        y(1:end)/delta,'--','color',plotSpecs.Cbudget(3,:),'linewidth',plotSpecs.linewidth)
    
    plot(dataWC(runNum).waveDiffusion(1:end,phase)/xnorm,...
        y(1:end)/delta,'color',plotSpecs.Cbudget(1,:),'linewidth',plotSpecs.linewidth)
    plot(dataWC(runNum).waveDeposition(1:end,phase)/xnorm,...
        y(1:end)/delta,'color',plotSpecs.Cbudget(2,:),'linewidth',plotSpecs.linewidth)
    plot(dataWC(runNum).waveSedTurb(1:end,phase)/xnorm,...
        y(1:end)/delta,'color',plotSpecs.Cbudget(3,:),'linewidth',plotSpecs.linewidth)
    
    set(gca,'FontSize',plotSpecs.axesfontsize)
    count = count+1;
end

% legend
clear h
h(1) = plot(0,0,'color','k');%, 'visible', 'off');
h(2) = plot(0,0,'--','color',plotSpecs.Cbudget(1,:));%, 'visible', 'off');
h(3) = plot(0,0,'--','color',plotSpecs.Cbudget(2,:));%, 'visible', 'off');
h(4) = plot(0,0,'--','color',plotSpecs.Cbudget(3,:));%, 'visible', 'off');
leg = legend(h,'$\frac{\partial \langle \widetilde{C} \rangle_{\hbox{p}}}{\partial t}$',...
    '$K \frac{\partial^2 C_{\hbox{c}}}{\partial z^2}$',...
    '$\frac{\partial}{\partial z} (w_{s} C_{\hbox{c}})$',...
    '-$\frac{\partial}{\partial z} (w''C'')_{\hbox{c}}$');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.836922805649892 0.860887412797838 0.120697937692915 0.105879838126046])

% Add text panel label for each panel
addText('$\theta = 0$', fig3,[0.308421502847779 0.928197022251018 0.0466913870402745 0.0295238091832115])
addText('$\theta = 1/6 \pi$', fig3,[0.549206752489743 0.922482736309666 0.100835173470633 0.0409523810659136])
addText('$\theta = 1/3 \pi$', fig3,[0.689456004809072 0.928197022251019 0.0667652402605329 0.0295238091832115])
addText('$\theta = 1/2 \pi$', fig3,[0.291241719094793 0.706768450822446 0.0667652402605329 0.0295238091832115])
addText('$\theta = 2/3 \pi$', fig3,[0.566241719094793 0.706768450822447 0.0667652402605329 0.0295238091832115])
addText('$\theta = 5/6 \pi$', fig3,[0.841241719094792 0.706768450822448 0.0667652402605329 0.0295238091832115])
addText('$\theta = \pi$', fig3,[0.309545319474192 0.485339879393875 0.048015182358878 0.0295238091832115])
addText('$\theta = 7/6 \pi$', fig3,[0.566241719094792 0.482958927012923 0.0667652402605329 0.0295238091832115])
addText('$\theta = 4/3 \pi$', fig3,[0.841241719094792 0.485339879393876 0.0667652402605329 0.0295238091832115])
addText('$\theta = 3/2 \pi$', fig3,[0.291241719094793 0.263911307965305 0.0667652402605329 0.0295238091832116])
addText('$\theta = 5/3 \pi$', fig3,[0.566241719094776 0.261530355584353 0.0667652402605329 0.0295238091832116])
addText('$\theta = 11/6 \pi$', fig3,[0.836330996294129 0.261530355584352 0.0730152572904315 0.0295238091832116])

if savePlot
    if runNum == 2
        print('./Figures/component_phaseSedBudget200WC_top','-depsc')
        print('./jpegs/component_phaseSedBudget200WC_top','-djpeg','-r600')
    elseif runNum == 4
        print('./Figures/component_phaseSedBudget350WC_top','-depsc')
        print('./jpegs/component_phaseSedBudget350WC_top','-djpeg','-r600')
    elseif runNum == 6
        print('./Figures/component_phaseSedBudgetWC_top','-depsc')
        print('./jpegs/component_phaseSedBudget500WC_top','-djpeg','-r600')
    end
end


%% Sediment fluxes integrated over bottom 5 delta in layers with thickness delta
Cmean = 1/H*dy*dataWC(6).CsedSteady; % depth-averaged sediment concentration used for flux normalization
deltaZ = 1*delta; % thickness of layers
fluxNorm = dataWC(6).ws*Cmean/deltaZ; % flux nomalization
integrationLim = [y(1)/delta, 1;1,2;2,3;3,4;4,5]; % integration limits in delta for each layer (1st column gives start height and 2nd column the end height
ylimAll = [-3, 7; -1, 1.5; -.3, .5; -.12, .2; -0.06, 0.06]; % y-axis limits for each layer

%%%%%%%Need to comment this function if I use it%%%%%%%%%%%
plotYYSedFluxOneRun(dataWC,y,delta,plotSpecs,6,integrationLim,ylimAll,0)

%plotYYSedFluxOneRun(dataWC,y,delta,plotSpecs,6,integrationLim,ylimAll,comp)

%legend
leg = legend('$\frac{\partial \langle \widetilde{C}_{\mathrm{layer}} \rangle_{\mathrm{p}}}{\partial t}$',...
    '$\Delta F_\mathrm{D}$','$\Delta F_\mathrm{w}$','$\Delta F_\mathrm{T}$');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.442438935382021 0.91859452469009 0.388237106800079 0.0380846619606018],...
    'Orientation','horizontal')
if savePlot
    print('./Figures/integratedFlux_500WC_bottom','-depsc')
    print('./jpegs/integratedFlux_500WC_bottom','-djpeg','-r600')
end

%% Panel plot of bed shear stress, fluxes, and sediment concentration
plotRun = 6; %Run to plot

% plots figure with five panels: 1) pressure and free-stream velocity
% 2) bed shear, 3) fluxes integrated from the bed to the top of the buffer
% layer, 4) fluxes integrated from the top of the buffer
% layer to the free-surface, 5) sediment conentration in each layer.
TwoLayerPanel(dataWC,y,dy,plotSpecs,params,plotRun)

if savePlot
    print('./Figures/TwoLayer_500WC','-depsc')
    print('./jpegs/TwoLayer_500WC','-djpeg','-r600')
end

%% Panel plot of bed shear stress and fluxes at the bed, top of viscous sublayer, and top of buffer layer
plotRun = 6; %Run to plot

% plots figure with 4-panels: 1) bed shear stress, 2-4) sediment fluxes at
% the bed, top of viscous sublayer, and buffer layer. Also output index for
% each plotting location.
[bedInd,viscousInd,bufferInd] = fixedZFluxes(dataWC,y,dy,plotSpecs,params,plotRun);

if savePlot
    print('./Figures/sedFluxFixedZ','-depsc')
    print('./jpegs/sedFluxFixedZ','-djpeg','-r600')
end

%% Panel plot of concentration at 3 hieghts were
% plots figure with 4-panels: 1) bed shear stress, 2-4) sediment fluxes at
% the bed, top of viscous sublayer, and buffer layer. Also output index for
% each plotting location.

fig8 = figure;
fig8.PaperUnits = 'centimeters';
fig8.PaperPosition = [0 0 14 8];
clf

xlim = [0 2*pi]; % x limits
ha = tight_subplot(3,1,[0.04 .025],[.12 .04],[.12 .08]);
panelCount = 1;
phasePlot = (0:num_steps_period)/num_steps_period*2*pi; % phases to plot over

for ind = [bedInd,viscousInd,bufferInd] % loop through height indexes to plot each
    colorCount =1; % counter for line color
    
    axes(ha(panelCount))
    set(gca,'box','on')
    hold
    for plotRun = [7, runsCup]
        phaseAveragedCsed =postProcess(dataWC(plotRun).Csed,0,y,0,'phaseAverage')*1000; % phase average the shear stress
        if plotRun == 7
           plot(phasePlot,[phaseAveragedCsed(ind,:), phaseAveragedCsed(ind,1)],'k','linewidth',plotSpecs.linewidth)         
        else
        plot(phasePlot,[phaseAveragedCsed(ind,:), phaseAveragedCsed(ind,1)],'color',plotSpecs.Cbudget(colorCount,:),'linewidth',plotSpecs.linewidth)
        colorCount =colorCount+1;
        end
    end
    if panelCount == 2 % middle panel
        set(gca,'xticklabel','','xtick', xtick,'xlim',xlim,'FontSize',plotSpecs.axesfontsize)
        myYlabel('$\langle C \rangle_{\mathrm{p}}$ (mg/L)',plotSpecs.labfontsize)
    elseif panelCount == 3 % bottom panel
        myXlabel('$\theta$',plotSpecs.labfontsize)
        set(gca,'xticklabel',xticklabels,'xtick', xtick,'xlim',xlim,'FontSize',plotSpecs.axesfontsize)
    elseif panelCount == 1 % top panel
        set(gca,'xticklabel','','xtick', xtick,'xlim',xlim,'FontSize',plotSpecs.axesfontsize)
        % legend
        leg = legend('0W','200WC', '350WC','500WC');
        set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
            'Position',[0.410020526090741 0.894238816909126 0.416716572216579 0.0307142853736877],...
            'Orientation','horizontal')
    end
    panelCount = panelCount+1;
end
addText('(a)', fig8,[0.125715569140678 0.913911307965305 0.0335318258830479 0.0295238091832115])
addText('(b)', fig8,[0.125368345688518 0.621054165108162 0.0342262727873666 0.0295238091832115])
addText('(c)', fig8,[0.124277087673567 0.328197022251022 0.0328373602458409 0.0295238091832116])


if savePlot
    print('./Figures/CsedFixedHeight','-depsc')
    print('./jpegs/CsedFixedHeight','-djpeg','-r600')
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

