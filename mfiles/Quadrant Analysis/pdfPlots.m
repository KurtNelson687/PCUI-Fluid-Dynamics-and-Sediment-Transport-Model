%%
close all; clear all; clear all
addpath('../Functions');

%% Variables
phase = 1;
plotHeight = 5;
runs = [1,2]; % Runs to plot
savePlot = true;
rhoWater = 1000;
thres = 3;
data_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/';
simulationType = 'strat200_1';
%% Data
load('../DataFolder/y') % loads y data
load('./figureProperties.mat') % loads figure properties for plotting
load('../DataFolder/CombinedHistData') % loads histogram data
load('../DataFolder/strat200.mat') % for uvRey

%% Find index for height
working_folder = [data_folder simulationType];


% -------------------------------------------------------------------------
% Get problem parameters and variables from PCUI
% -------------------------------------------------------------------------
% read the file containing the parameter definition
ftext = fileread(fullfile(working_folder, 'pcuiRunParams.txt'));
params.molecular_viscosity = variable_value_pcui('vis',ftext);
params.dpdxSteady = variable_value_pcui('dpdxSteady',ftext);
params.rhoSed = variable_value_pcui('rhoSed',ftext);
H = y(end)+(y(end)-y(end-1))/2;
ustar = sqrt(params.dpdxSteady*H/rhoWater);
UVindex = find(y/(params.molecular_viscosity/ustar) >= plotHeight,1,'first');


%% Plot pdfs
fig1 = figure;
fig1.PaperUnits = 'centimeters';
fig1.PaperPosition = [0 0 14 7];
clf
ha = tight_subplot(1,2,[0.03 .04],[.12 .06],[.08 .06]);
count = 1;
ylim = [-6*10^-4 6*10^-4];
xlim = [-0.01 0.01];

plotCount = 1;
for plotNum = runs
    heightInd = find(allHeightData(plotNum).pdfHeights==plotHeight); % index for plot height
    
    % plotting
    axes(ha(plotCount))
    set(gca,'box','on')
    his_plot = allHistData(plotNum,heightInd,phase).uvMatrix; % extract pdf
    numTot = sum(sum(his_plot)); % total number of events
    his_plot(his_plot == 0) = NaN;
    his_plot = his_plot/numTot;
    u = allEdgeData(plotNum,heightInd).uEdge;
    v = allEdgeData(plotNum,heightInd).vEdge;
    h = pcolor(u(1:end-1),v(1:end-1),his_plot);
    set(h, 'EdgeColor', 'none');
    set(gca,'FontSize',plotSpecs.axesfontsize,'ylim',ylim,'xlim',xlim)
    
    
    if plotCount == 1
        check = nansum(nansum(his_plot))
        cmax = max(max(his_plot));
        ylab = ylabel('$w \prime$');
        xlab = xlabel('$u \prime$');
        set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
        set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
        %addText('(a)', fig1,[0.139091229915617 0.905476190646492 0.0335318258830479 0.0295238091832115])
    elseif plotCount == 2
        check = nansum(nansum(his_plot))
       % addText('(b)', fig1,[0.548019801344188 0.905476190646493 0.0335318258830479 0.0295238091832115])
        set(gca,'yticklabel','')
        xlab = xlabel('$u \prime$');
        set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
    end
    
    hold
    plotThres(xlim,ylim,thres,uvReyPhaseAve(UVindex))

    caxis([0 cmax])
    plotCount = plotCount +1;
end



if savePlot
    print('./Figures/uvpdf','-depsc')
    print('./jpegs/uvpdf','-djpeg','-r600')
end

%% Plot uvpdf diff
fig2 = figure;
fig2.PaperUnits = 'centimeters';
fig2.PaperPosition = [0 0 9 7];
set(gca,'box','on')
%ylim = [-3*10^-3 3*10^-3];
%xlim = [-0.03 0.03];

heightInd = find(allHeightData(plotNum).pdfHeights==plotHeight); % index for plot height

set(gca,'box','on')
his_noCup = allHistData(runs(1),heightInd,phase).uvMatrix; % extract pdf
numTot_nocup = sum(sum(his_noCup)); % total number of events
his_noCup = his_noCup/numTot_nocup;

his_Cup = allHistData(runs(2),heightInd,phase).uvMatrix; % extract pdf
numTot_cup = sum(sum(his_Cup)); % total number of events
his_Cup = his_Cup/numTot_cup;

his_plot = his_noCup-his_Cup;
his_plot(his_plot == 0) = NaN;
u = allEdgeData(runs(1),heightInd).uEdge;
v = allEdgeData(runs(1),heightInd).vEdge;
h = pcolor(u(1:end-1),v(1:end-1),his_plot);
set(h, 'EdgeColor', 'none');
set(gca,'FontSize',plotSpecs.axesfontsize,'ylim',ylim,'xlim',xlim)
ylab = ylabel('$w \prime$');
xlab = xlabel('$u \prime$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
cmax = max(max(his_plot));
cmin = min(min(his_plot));
caxis([cmin cmax])
colorbar

hold
plotThres(xlim,ylim,thres,uvReyPhaseAve(UVindex))




if savePlot
    print('./Figures/uvpdfDiff','-depsc')
    print('./jpegs/uvpdfDiff','-djpeg','-r600')
end

%% Plot pdfs who
fig3 = figure;
fig3.PaperUnits = 'centimeters';
fig3.PaperPosition = [0 0 14 7];
clf
ha = tight_subplot(1,2,[0.03 .04],[.12 .06],[.08 .06]);
count = 1;
ylim = [-6*10^-4 6*10^-4];
xlim = [-0.01 0.01];

plotCount = 1;
for plotNum = runs
    heightInd = find(allHeightData(plotNum).pdfHeights==plotHeight); % index for plot height
    
    % plotting
    axes(ha(plotCount))
    set(gca,'box','on')
    his_plot = allHistData(plotNum,heightInd,phase).vCsedMatrix; % extract pdf
    numTot = sum(sum(his_plot)); % total number of events
    his_plot(his_plot == 0) = NaN;
    his_plot = his_plot/numTot;
    u = allEdgeData(plotNum,heightInd).rhoEdge; 
    v = allEdgeData(plotNum,heightInd).vEdge;
    h = pcolor(u(1:end-1),v(1:end-1),his_plot);
    set(h, 'EdgeColor', 'none');
    set(gca,'FontSize',plotSpecs.axesfontsize,'ylim',ylim)%,'xlim',xlim)
    
    
    if plotCount == 1
        check = nansum(nansum(his_plot))
        cmax = max(max(his_plot));
        ylab = ylabel('$w \prime$');
        xlab = xlabel('$\rho \prime$');
        set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
        set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
        %addText('(a)', fig1,[0.139091229915617 0.905476190646492 0.0335318258830479 0.0295238091832115])
    elseif plotCount == 2
        check = nansum(nansum(his_plot))
        %addText('(b)', fig1,[0.548019801344188 0.905476190646493 0.0335318258830479 0.0295238091832115])
        set(gca,'yticklabel','')
        xlab = xlabel('$\rho \prime$');
        set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
    end
    
    hold
    plotThres(xlim,ylim,thres,uvReyPhaseAve(UVindex))

    caxis([0 cmax])
    plotCount = plotCount +1;
end



if savePlot
    print('./Figures/vrhopdf','-depsc')
    print('./jpegs/vrhopdf','-djpeg','-r600')
end

%% Plot uvpdf diff
fig4 = figure;
fig4.PaperUnits = 'centimeters';
fig4.PaperPosition = [0 0 9 7];
set(gca,'box','on')
%ylim = [-3*10^-3 3*10^-3];
%xlim = [-0.03 0.03];

heightInd = find(allHeightData(plotNum).pdfHeights==plotHeight); % index for plot height

set(gca,'box','on')
his_noCup = allHistData(runs(1),heightInd,phase).vCsedMatrix; % extract pdf
numTot_nocup = sum(sum(his_noCup)); % total number of events
his_noCup = his_noCup/numTot_nocup;

his_Cup = allHistData(runs(2),heightInd,phase).vCsedMatrix; % extract pdf
numTot_cup = sum(sum(his_Cup)); % total number of events
his_Cup = his_Cup/numTot_cup;

his_plot = his_noCup-his_Cup;
his_plot(his_plot == 0) = NaN;
u = allEdgeData(runs(1),heightInd).rhoEdge;
v = allEdgeData(runs(1),heightInd).vEdge;
h = pcolor(u(1:end-1),v(1:end-1),his_plot);
set(h, 'EdgeColor', 'none');
set(gca,'FontSize',plotSpecs.axesfontsize,'ylim',ylim,'xlim',xlim)
ylab = ylabel('$w \prime$');
xlab = xlabel('$\rho \prime$');
set(ylab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
set(xlab,'interpreter','Latex','FontSize',plotSpecs.labfontsize)
cmax = max(max(his_plot));
cmin = min(min(his_plot));
caxis([cmin cmax])
colorbar

hold
plotThres(xlim,ylim,thres,uvReyPhaseAve(UVindex))




if savePlot
    print('./Figures/vrhopdfDiff','-depsc')
    print('./jpegs/vrhopdfDiff','-djpeg','-r600')
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

function plotThres(xlim,ylim,thres,uv)
thresVal = thres*uv;
u = linspace(xlim(1),xlim(2),300);
v1 = thresVal./u;
v2 = -1*thresVal./u;

plot(xlim,[0, 0],'k')
plot([0, 0],ylim,'k')
plot(u,v1,':k')
plot(u,v2,':k')
set(gca,'ylim',ylim,'xlim',xlim)
end

