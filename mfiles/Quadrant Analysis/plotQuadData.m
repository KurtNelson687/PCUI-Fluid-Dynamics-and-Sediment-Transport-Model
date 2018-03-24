%%
close all; clear all; clear all
addpath('../Functions');

%% Variables
runs = [1,2]; % Runs to plot
% x labels
xticklabels = {'0','\pi/4', '\pi/2','3\pi/4', '\pi','5\pi/4', '3\pi/2','7\pi/4', '2 \pi'};
xtick = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4,2*pi];
savePlot = true;
%% Data
load('../DataFolder/CombinedQuadData')
load('../DataFolder/y')
load('./figureProperties.mat') %figure properties for plotting

%% Plot Reynolds stress contribution with height - period averaged
fig1 = figure;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0 0 14 15];
clf
ha = tight_subplot(3,1,[0.03 .02],[.1 .04],[.12 .08]);
count = 1;


% panel 1
ylim = [-2 3.5];
xlim = [0 30];
axes(ha(1))
set(gca,'box','on')
hold
thres = 0;
thresInd = find(allThresholds>=thres,1,'first');
for plotNum = runs
    plotFixThresPeriodQ(allQuadData,allHeightData(plotNum).allHeights,plotNum,plotSpecs,thresInd)
end
set(gca,'FontSize',plotSpecs.axesfontsize,'ylim',ylim,'xlim',xlim)
set(gca,'xticklabel','')

leg = legend('Q1','Q2','Q3','Q4');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.464212049408042 0.923180428889602 0.340674795423235 0.0307142853736877],...
    'Orientation','horizontal')

addText('Total', fig1,[0.731357142857141 0.707142857142858 0.141857142857143 0.0309523809523817])

% panel 2
ylim = [-1 1.5];
axes(ha(2))
set(gca,'box','on')
hold
thres = 3;
thresInd = find(allThresholds>=thres,1,'first');
for plotNum = runs
    plotFixThresPeriodQSplit(allQuadData,allHeightData(plotNum).allHeights,plotNum,plotSpecs,thresInd,1)
end
set(gca,'FontSize',plotSpecs.axesfontsize,'ylim',ylim,'xlim',xlim)
set(gca,'xticklabel','')

addText('Below threshold', fig1,[0.716071428571429 0.407142857142858 0.212499999999998 0.0357142857142856])

% panel 3
axes(ha(3))
set(gca,'box','on')
hold
for plotNum = runs
    plotFixThresPeriodQSplit(allQuadData,allHeightData(plotNum).allHeights,plotNum,plotSpecs,thresInd,2)
end
set(gca,'FontSize',plotSpecs.axesfontsize,'ylim',ylim,'xlim',xlim)
xlab = xlabel('$z^{+}$');
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)

addText('Above threshold', fig1,[0.722428571428568 0.100000000000004 0.209714285714287 0.0428571428571595])

if savePlot
    print('./Figures/stressContributionProfiles','-depsc')
    print('./jpegs/stressContributionProfiles','-djpeg','-r600')
end

%% Plot Reynolds stress contribution with height - period averaged
fig4 = figure;
fig4.PaperUnits = 'centimeters';
fig4.PaperPosition = plotSpecs.fullWidthHalfPage;
clf
ha = tight_subplot(2,1,[0.03 .02],[.1 .04],[.12 .08]);
count = 1;
ylim = [0 1];

% panel 1
axes(ha(1))
set(gca,'box','on')
hold
Hplot = 20;
for plotNum = runs
    Hindex = find(allHeightData(plotNum).allHeights>=Hplot,1,'first');
    plotFixHPeriodN(allQuadData,Hindex,plotNum,plotSpecs,allThresholds)
end
set(gca,'FontSize',plotSpecs.axesfontsize)
set(gca,'xticklabel','')

leg = legend('Q1','Q2','Q3','Q4');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.464212049408042 0.923180428889602 0.340674795423235 0.0307142853736877],...
    'Orientation','horizontal')

% panel 2
axes(ha(2))
set(gca,'box','on')
hold
Hplot = 30;
for plotNum = runs
    Hindex = find(allHeightData(plotNum).allHeights>=Hplot,1,'first');
    plotFixHPeriodN(allQuadData,Hindex,plotNum,plotSpecs,allThresholds)
end
set(gca,'FontSize',plotSpecs.axesfontsize)
xlab = xlabel('Threshold');
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)

if savePlot
    print('./jpegs/StressQuadrantCount_FixedH','-djpeg','-r600')
end

%% Plot Reynolds stress contribution fixed height - period averaged
fig5 = figure;
fig5.PaperUnits = 'centimeters';
fig5.PaperPosition = plotSpecs.fullWidthHalfPage;
clf
ha = tight_subplot(2,1,[0.03 .02],[.1 .04],[.12 .08]);
count = 1;
ylim = [-2.5 4];
xlim = [0 2*pi];

% panel 1
axes(ha(1))
set(gca,'box','on')
hold
Hplot = 30;
thres = 5;
for plotNum = runs
    thresInd = find(allThresholds>=thres,1,'first');
    Hindex = find(allHeightData(plotNum).allHeights>=Hplot,1,'first');
    plotFixHPhaseQ(allQuadData,Hindex,plotNum,plotSpecs,thresInd)
end
set(gca,'FontSize',plotSpecs.axesfontsize,'ylim',ylim,'xlim',xlim)
set(gca,'xtick',xtick,'xticklabel','')

ylab = ylabel('Stress Fraction');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)


leg = legend('Q1','Q2','Q3','Q4');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.464212049408042 0.923180428889602 0.340674795423235 0.0307142853736877],...
    'Orientation','horizontal')

% panel 2
ylim = [0 1];
axes(ha(2))
set(gca,'box','on')
hold
Hplot = 20;
thres = 5;
for plotNum = runs
    thresInd = find(allThresholds>=thres,1,'first');
    Hindex = find(allHeightData(plotNum).allHeights>=Hplot,1,'first');
    plotFixHPhaseN(allQuadData,Hindex,plotNum,plotSpecs,thresInd)
end
ylab = ylabel('Counts');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)

set(gca,'FontSize',plotSpecs.axesfontsize,'xlim',xlim,'xtick',xtick,'xticklabels',xticklabels)
xlab = xlabel('Threshold');
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)

if savePlot
    print('./jpegs/stressPhase','-djpeg','-r600')
end

%% Functions
function plotFixThresPeriodQ(allQuadData,allHeights,plotNum,plotSpecs,thresInd)
if plotNum == 1
    load('../dataForQuad/strat200.mat')
    working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/strat200_2';
elseif plotNum == 2
    load('../dataForQuad/strat200cup.mat')
    working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/strat200_2';
elseif plotNum == 3
    load('../dataForQuad/strat350.mat')
    working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/strat350_2';
elseif plotNum == 4
    load('../dataForQuad/strat350cup.mat')
    working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/strat350_2';
elseif plotNum == 5
    load('../dataForQuad/strat500.mat')
    working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/strat500_2';
elseif plotNum == 6
    load('../dataForQuad/strat500cup.mat')
    working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/strat500_2';
end

%%%%%%%%%%%%%%%%%%%%%%%%% For QC %%%%%%%%%%%%%%%
rhoWater = 1000;
load('../DataFolder/y')
H = y(end)+(y(end)-y(end-1))/2;
ftext = fileread(fullfile(working_folder, 'pcuiRunParams.txt'));
params.molecular_viscosity = variable_value_pcui('vis',ftext);
params.dpdxSteady = variable_value_pcui('dpdxSteady',ftext);
params.rhoSed = variable_value_pcui('rhoSed',ftext);
ustar = sqrt(params.dpdxSteady*H/rhoWater);
allHeights = [1:60,70:20:max(y/(params.molecular_viscosity/ustar))];

ustar = sqrt(params.dpdxSteady*H/rhoWater);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(allHeights)
    Hindex_temp = find(y/(params.molecular_viscosity/ustar) >= allHeights(i),1,'first');
    
    if mod(plotNum,2) == 0
        normNum = plotNum;
        normThres = thresInd;
    else
        normNum = plotNum;
        normThres = thresInd;
    end
    
    numPoints = sum(allQuadData(plotNum,i,normThres).Nstress(1,:)) + ...
        sum(allQuadData(plotNum,i,normThres).Nstress(2,:)) + ...
        sum(allQuadData(plotNum,i,normThres).Nstress(3,:)) + ...
        sum(allQuadData(plotNum,i,normThres).Nstress(4,:));
    
    numPointsNorm = sum(allQuadData(normNum,i,normThres).Nstress(1,:)) + ...
        sum(allQuadData(normNum,i,normThres).Nstress(2,:)) + ...
        sum(allQuadData(normNum,i,normThres).Nstress(3,:)) + ...
        sum(allQuadData(normNum,i,normThres).Nstress(4,:));
    
    tot =   (sum(allQuadData(normNum,i,normThres).Qstress(1,:)) +...
        sum(allQuadData(normNum,i,normThres).Qstress(2,:)) +...
        sum(allQuadData(normNum,i,normThres).Qstress(3,:)) +...
        sum(allQuadData(normNum,i,normThres).Qstress(4,:)))./sum(numPointsNorm);
    
    ratio = tot/uvReyPhaseAve(Hindex_temp);
    
    if abs(ratio)<0.8
        'small ratio'
        'runNum'
        plotNum
        'height'
        allHeights(i)
    end
    
    plotVar1(i) = sum(allQuadData(plotNum,i,thresInd).Qstress(1,:))/sum(numPoints)/tot;
    plotVar2(i) = sum(allQuadData(plotNum,i,thresInd).Qstress(2,:))/sum(numPoints)/tot;
    plotVar3(i) = sum(allQuadData(plotNum,i,thresInd).Qstress(3,:))/sum(numPoints)/tot;
    plotVar4(i) = sum(allQuadData(plotNum,i,thresInd).Qstress(4,:))/sum(numPoints)/tot;
    
end

if mod(plotNum,2) == 0
    linetype = '--';
else
    linetype = '-';
end

plotStart = 1;
semilogx(allHeights(plotStart:end),plotVar1(plotStart:end),linetype,'color','k','linewidth',plotSpecs.linewidth)
semilogx(allHeights(plotStart:end),plotVar2(plotStart:end),linetype,'color','r','linewidth',plotSpecs.linewidth)
semilogx(allHeights(plotStart:end),plotVar3(plotStart:end),linetype,'color','b','linewidth',plotSpecs.linewidth)
semilogx(allHeights(plotStart:end),plotVar4(plotStart:end),linetype,'color','g','linewidth',plotSpecs.linewidth)
end

function plotFixThresPeriodN(allQuadData,allHeights,plotNum,plotSpecs,thresInd)

for i = 1:length(allHeights)
    if mod(plotNum,2) == 0
        normNum = plotNum;
    else
        normNum = plotNum;
    end
    
    tot =   sum(allQuadData(normNum,i,1).Nstress(1,:)) +...
        sum(allQuadData(normNum,i,1).Nstress(2,:)) +...
        sum(allQuadData(normNum,i,1).Nstress(3,:)) +...
        sum(allQuadData(normNum,i,1).Nstress(4,:));
    
    plotVar1(i) = sum(allQuadData(plotNum,i,thresInd).Nstress(1,:))/tot;
    plotVar2(i) = sum(allQuadData(plotNum,i,thresInd).Nstress(2,:))/tot;
    plotVar3(i) = sum(allQuadData(plotNum,i,thresInd).Nstress(3,:))/tot;
    plotVar4(i) = sum(allQuadData(plotNum,i,thresInd).Nstress(4,:))/tot;
end

if mod(plotNum,2) == 0
    linetype = '--';
else
    linetype = '-';
end

semilogx(allHeights,plotVar1,linetype,'color','k','linewidth',plotSpecs.linewidth)
semilogx(allHeights,plotVar2,linetype,'color','r','linewidth',plotSpecs.linewidth)
semilogx(allHeights,plotVar3,linetype,'color','b','linewidth',plotSpecs.linewidth)
semilogx(allHeights,plotVar4,linetype,'color','g','linewidth',plotSpecs.linewidth)
end

function plotFixHPeriodQ(allQuadData,Hindex,plotNum,plotSpecs,allThresholds)
for i = 1:length(allThresholds)
    tot =   sum(allQuadData(plotNum,Hindex,i).Qstress(1,:)) + ...
        sum(allQuadData(plotNum,Hindex,i).Qstress(2,:)) + ...
        sum(allQuadData(plotNum,Hindex,i).Qstress(3,:)) + ...
        sum(allQuadData(plotNum,Hindex,i).Qstress(4,:))./...
        (sum(allQuadData(plotNum,Hindex,i).Nstress(1,:)) + ...
        sum(allQuadData(plotNum,Hindex,i).Nstress(2,:)) + ...
        sum(allQuadData(plotNum,Hindex,i).Nstress(3,:)) + ...
        sum(allQuadData(plotNum,Hindex,i).Nstress(4,:)));
    
    
    plotVar1(i) = sum(allQuadData(plotNum,Hindex,i).Qstress(1,:))./...
        sum(allQuadData(plotNum,Hindex,i).Nstress(1,:))./tot;
    plotVar2(i) = sum(allQuadData(plotNum,Hindex,i).Qstress(2,:))./...
        sum(allQuadData(plotNum,Hindex,i).Nstress(2,:))./tot;
    plotVar3(i) = sum(allQuadData(plotNum,Hindex,i).Qstress(3,:))./...
        sum(allQuadData(plotNum,Hindex,i).Nstress(3,:))./tot;
    plotVar4(i) = sum(allQuadData(plotNum,Hindex,i).Qstress(4,:))./...
        sum(allQuadData(plotNum,Hindex,i).Nstress(4,:))./tot;
end

if mod(plotNum,2) == 0
    linetype = '--';
else
    linetype = '-';
end

plot(allThresholds,plotVar1,linetype,'color','k','linewidth',plotSpecs.linewidth)
plot(allThresholds,plotVar2,linetype,'color','r','linewidth',plotSpecs.linewidth)
plot(allThresholds,plotVar3,linetype,'color','b','linewidth',plotSpecs.linewidth)
plot(allThresholds,plotVar4,linetype,'color','g','linewidth',plotSpecs.linewidth)
end

function plotFixHPeriodN(allQuadData,Hindex,plotNum,plotSpecs,allThresholds)
for i = 1:length(allThresholds)
    if mod(plotNum,2) == 0
        normNum = plotNum;
    else
        normNum = plotNum;
    end
    
    tot =   sum(allQuadData(normNum,Hindex,1).Nstress(1,:)) + ...
        sum(allQuadData(normNum,Hindex,1).Nstress(2,:)) + ...
        sum(allQuadData(normNum,Hindex,1).Nstress(3,:)) + ...
        sum(allQuadData(normNum,Hindex,1).Nstress(4,:));
    
    
    plotVar1(i) = sum(allQuadData(plotNum,Hindex,i).Nstress(1,:))/tot;
    plotVar2(i) = sum(allQuadData(plotNum,Hindex,i).Nstress(2,:))/tot;
    plotVar3(i) = sum(allQuadData(plotNum,Hindex,i).Nstress(3,:))/tot;
    plotVar4(i) = sum(allQuadData(plotNum,Hindex,i).Nstress(4,:))/tot;
end

if mod(plotNum,2) == 0
    linetype = '--';
else
    linetype = '-';
end

plot(allThresholds,plotVar1,linetype,'color','k','linewidth',plotSpecs.linewidth)
plot(allThresholds,plotVar2,linetype,'color','r','linewidth',plotSpecs.linewidth)
plot(allThresholds,plotVar3,linetype,'color','b','linewidth',plotSpecs.linewidth)
plot(allThresholds,plotVar4,linetype,'color','g','linewidth',plotSpecs.linewidth)
end

function plotFixHPhaseQ(allQuadData,Hindex,plotNum,plotSpecs,thresInd)
xvar = 0:pi/6:2*pi;

numPoints = (allQuadData(plotNum,Hindex,thresInd).Nstress(1,:)) + ...
    sum(allQuadData(plotNum,Hindex,thresInd).Nstress(2,:)) + ...
    sum(allQuadData(plotNum,Hindex,thresInd).Nstress(3,:)) + ...
    sum(allQuadData(plotNum,Hindex,thresInd).Nstress(4,:));

tot =   (allQuadData(plotNum,Hindex,thresInd).Qstress(1,:) + ...
    allQuadData(plotNum,Hindex,thresInd).Qstress(2,:) + ...
    allQuadData(plotNum,Hindex,thresInd).Qstress(3,:) + ...
    allQuadData(plotNum,Hindex,thresInd).Qstress(4,:))/numPoints;


plotVar1 = allQuadData(plotNum,Hindex,thresInd).Qstress(1,:)./numPoints./tot;
plotVar2 = allQuadData(plotNum,Hindex,thresInd).Qstress(2,:)./numPoints./tot;
plotVar3 = allQuadData(plotNum,Hindex,thresInd).Qstress(3,:)./numPoints./tot;
plotVar4 = allQuadData(plotNum,Hindex,thresInd).Qstress(4,:)./numPoints./tot;

if mod(plotNum,2) == 0
    linetype = '--';
else
    linetype = '-';
end

plot(xvar,[plotVar1, plotVar1(1)],linetype,'color','k','linewidth',plotSpecs.linewidth)
plot(xvar,[plotVar2, plotVar2(1)],linetype,'color','r','linewidth',plotSpecs.linewidth)
plot(xvar,[plotVar3, plotVar3(1)],linetype,'color','b','linewidth',plotSpecs.linewidth)
plot(xvar,[plotVar4, plotVar4(1)],linetype,'color','g','linewidth',plotSpecs.linewidth)
end

function plotFixHPhaseN(allQuadData,Hindex,plotNum,plotSpecs,thresInd)
xvar = 0:pi/6:2*pi;

%tot = 1;
tot =   allQuadData(plotNum,Hindex,thresInd).Nstress(1,:) + ...
    allQuadData(plotNum,Hindex,thresInd).Nstress(2,:) + ...
    allQuadData(plotNum,Hindex,thresInd).Nstress(3,:) + ...
    allQuadData(plotNum,Hindex,thresInd).Nstress(4,:);


plotVar1 = allQuadData(plotNum,Hindex,thresInd).Nstress(1,:)./tot;
plotVar2 = allQuadData(plotNum,Hindex,thresInd).Nstress(2,:)./tot;
plotVar3 = allQuadData(plotNum,Hindex,thresInd).Nstress(3,:)./tot;
plotVar4 = allQuadData(plotNum,Hindex,thresInd).Nstress(4,:)./tot;

if mod(plotNum,2) == 0
    linetype = '--';
else
    linetype = '-';
end

plot(xvar,[plotVar1, plotVar1(1)],linetype,'color','k','linewidth',plotSpecs.linewidth)
plot(xvar,[plotVar2, plotVar2(1)],linetype,'color','r','linewidth',plotSpecs.linewidth)
plot(xvar,[plotVar3, plotVar3(1)],linetype,'color','b','linewidth',plotSpecs.linewidth)
plot(xvar,[plotVar4, plotVar4(1)],linetype,'color','g','linewidth',plotSpecs.linewidth)
end


function plotFixThresPeriodQSplit(allQuadData,allHeights,plotNum,plotSpecs,thresInd,direction)
if plotNum == 1
    load('../dataForQuad/strat200.mat')
    working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/strat200_2';
elseif plotNum == 2
    load('../dataForQuad/strat200cup.mat')
    working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/strat200_2';
elseif plotNum == 3
    load('../dataForQuad/strat350.mat')
    working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/strat350_2';
elseif plotNum == 4
    load('../dataForQuad/strat350cup.mat')
    working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/strat350_2';
elseif plotNum == 5
    load('../dataForQuad/strat500.mat')
    working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/strat500_2';
elseif plotNum == 6
    load('../dataForQuad/strat500cup.mat')
    working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/strat500_2';
end

%%%%%%%%%%%%%%%%%%%%%%%%% For QC %%%%%%%%%%%%%%%
rhoWater = 1000;
load('../DataFolder/y')
H = y(end)+(y(end)-y(end-1))/2;
ftext = fileread(fullfile(working_folder, 'pcuiRunParams.txt'));
params.molecular_viscosity = variable_value_pcui('vis',ftext);
params.dpdxSteady = variable_value_pcui('dpdxSteady',ftext);
params.rhoSed = variable_value_pcui('rhoSed',ftext);
ustar = sqrt(params.dpdxSteady*H/rhoWater);
allHeights = [1:60,70:20:max(y/(params.molecular_viscosity/ustar))];

ustar = sqrt(params.dpdxSteady*H/rhoWater);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(allHeights)
    Hindex_temp = find(y/(params.molecular_viscosity/ustar) >= allHeights(i),1,'first');
    
    if mod(plotNum,2) == 0
        normNum = plotNum;
    else
        normNum = plotNum;
    end
    
    
    numPointsNorm = sum(allQuadData(normNum,i,1).Nstress(1,:)) + ...
        sum(allQuadData(normNum,i,1).Nstress(2,:)) + ...
        sum(allQuadData(normNum,i,1).Nstress(3,:)) + ...
        sum(allQuadData(normNum,i,1).Nstress(4,:));
    
    tot =   (sum(allQuadData(normNum,i,1).Qstress(1,:)) +...
        sum(allQuadData(normNum,i,1).Qstress(2,:)) +...
        sum(allQuadData(normNum,i,1).Qstress(3,:)) +...
        sum(allQuadData(normNum,i,1).Qstress(4,:)))./sum(numPointsNorm);
    
    ratio = tot/uvReyPhaseAve(Hindex_temp);
    
    if abs(ratio)<0.8
        'small ratio'
        'runNum'
        plotNum
        'height'
        allHeights(i)
    end
    
    Above1(i) = sum(allQuadData(plotNum,i,thresInd).Qstress(1,:))./sum(numPointsNorm);
    Above2(i) = sum(allQuadData(plotNum,i,thresInd).Qstress(2,:))./sum(numPointsNorm);
    Above3(i) = sum(allQuadData(plotNum,i,thresInd).Qstress(3,:))./sum(numPointsNorm);
    Above4(i) = sum(allQuadData(plotNum,i,thresInd).Qstress(4,:))./sum(numPointsNorm);
    
    Below1(i) = sum(allQuadData(normNum,i,1).Qstress(1,:))./sum(numPointsNorm)-Above1(i);
    Below2(i) = sum(allQuadData(normNum,i,1).Qstress(2,:))./sum(numPointsNorm)-Above2(i);
    Below3(i) = sum(allQuadData(normNum,i,1).Qstress(3,:))./sum(numPointsNorm)-Above3(i);
    Below4(i) = sum(allQuadData(normNum,i,1).Qstress(4,:))./sum(numPointsNorm)-Above4(i);
    
    if direction == 1
        plotVar1(i) = Below1(i)/tot;
        plotVar2(i) = Below2(i)/tot;
        plotVar3(i) = Below3(i)/tot;
        plotVar4(i) = Below4(i)/tot;
    elseif direction == 2
        plotVar1(i) = Above1(i)/tot;
        plotVar2(i) = Above2(i)/tot;
        plotVar3(i) = Above3(i)/tot;
        plotVar4(i) = Above4(i)/tot;
    end
    
    
end

if mod(plotNum,2) == 0
    linetype = '--';
else
    linetype = '-';
end

plotStart = 1;
plot(allHeights(plotStart:end),plotVar1(plotStart:end),linetype,'color','k','linewidth',plotSpecs.linewidth)
plot(allHeights(plotStart:end),plotVar2(plotStart:end),linetype,'color','r','linewidth',plotSpecs.linewidth)
plot(allHeights(plotStart:end),plotVar3(plotStart:end),linetype,'color','b','linewidth',plotSpecs.linewidth)
plot(allHeights(plotStart:end),plotVar4(plotStart:end),linetype,'color','g','linewidth',plotSpecs.linewidth)
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


