%% Flux into top of water column
function [h] = plotYYSedFluxOneRun(dataWC,y,delta,plotSpecs,plotNum,integrationLim,ylimAll,comp)
% Function: plotYYSedFluxOneRun
%
% Author: Kurt Nelson
%
% Purpose: This function creates a 5-panel plot where each panel shows
% sediment fluxes integrated over a layer.
%
% Inputs:
% 1) dataWC - PCUI data array created by data extractor that contains
%    velocity and sediment information
% 2) y - vector containing cell y coordinate
% 3) delta - stokes layer thickness
% 4) plotSpecs - array containing plot specifications
% 5) plotNum - run to plot
% 6) integrationLim - integration limits for each layer
% 7) ylimAll - y limits for all plots
%
%%
fig = figure;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = plotSpecs.fullPageFigSize;
clf
ha = tight_subplot(5,1,[0.03 .025],[.1 .04],[.12 .08]);

% x labels
xticklabels = {'0','\pi/4', '\pi/2','3\pi/4', '\pi','5\pi/4', '3\pi/2','7\pi/4', '2 \pi'};
xtick = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4,2*pi];

fluxNorm = 10^-4; % flux normalization

for panel = 1:5
    set(gca,'box','on')
    
    if panel == 1 %panel 1
        axes(ha(1))
        bottomY = integrationLim(1,1); % bottom limit for integration for panel
        topY = integrationLim(1,2);  % top limit for integration for panel
        ylim = ylimAll(1,:); % y-limits for plot
        xlim = [0 2*pi]; % x-limits for plot
        set(gca,'ylim',ylim,'xlim',xlim,'xticklabel',[],'xtick',xtick,'FontSize',plotSpecs.axesfontsize) %axes properties
        %add text giving layer start and stop
        addText(['$z_1 =$ ' num2str(round(integrationLim(1,1))) '$\Delta$ to $z_1 =$ ' num2str(integrationLim(1,2)) '$\Delta$'],...
            fig, [0.766144283885229 0.892315012479271 0.126685796465193 0.0295238091832115]);
    elseif panel == 2 %panel 2 - same as panel 1 but different integration limits
        axes(ha(2))
        bottomY = integrationLim(2,1);
        topY = integrationLim(2,2);
        ylim = ylimAll(2,:);
        xlim = [0 2*pi];
        set(gca,'ylim',ylim,'xlim',xlim,'xticklabel',[],'xtick',xtick,'FontSize',plotSpecs.axesfontsize) %axes properties
        addText(['$z_1 =$ ' num2str(integrationLim(2,1)) '$\Delta$ to $z_1 =$ ' num2str(integrationLim(2,2)) '$\Delta$'],...
            fig, [0.764358569599515 0.73517215533642 0.126685796465193 0.0295238091832115]);
    elseif panel == 3 %panel 3 - same as panel 1 but different integration limits
        axes(ha(3))
        bottomY = integrationLim(3,1);
        topY = integrationLim(3,2);
        ylim = ylimAll(3,:);
        xlim = [0 2*pi];
        myYlabel('Sediment flux normalized by $w_{\mathrm{s}} \overline{C}/\Delta Z$',plotSpecs.labfontsize)
        set(gca,'ylim',ylim,'xlim',xlim,'xticklabel',[],'xtick',xtick,'FontSize',plotSpecs.axesfontsize) %axes properties
        addText(['$z_1 =$ ' num2str(integrationLim(3,1)) '$\Delta$ to $z_1 =$ ' num2str(integrationLim(3,2)) '$\Delta$'],...
            fig, [0.764358569599515 0.556600726764994 0.126685796465193 0.0295238091832115]);
    elseif panel == 4 %panel 4 - same as panel 1 but different integration limits
        axes(ha(4))
        bottomY = integrationLim(4,1);
        topY = integrationLim(4,2);
        ylim = ylimAll(4,:);
        xlim = [0 2*pi];
        set(gca,'ylim',ylim,'xlim',xlim,'xticklabel',[],'xtick',xtick,'FontSize',plotSpecs.axesfontsize) %axes properties
        addText(['$z_1 =$ ' num2str(integrationLim(4,1)) '$\Delta$ to $z_1 =$ ' num2str(integrationLim(4,2)) '$\Delta$'],...
            fig, [0.764358569599515 0.378029298193566 0.126685796465193 0.0295238091832116]);
    elseif panel ==5 %panel 5 - same as panel 1 but different integration limits
        axes(ha(5))
        bottomY = integrationLim(5,1);
        topY = integrationLim(5,2);
        ylim = ylimAll(5,:);
        xlim = [0 2*pi];
        set(gca,'ylim',ylim,'xlim',xlim,'xtick',xtick,'xticklabel',xticklabels,'FontSize',plotSpecs.axesfontsize)
        myXlabel('$\theta$',plotSpecs.labfontsize)
        addText(['$z_1 =$ ' num2str(integrationLim(5,1)) '$\Delta$ to $z_1 =$ ' num2str(integrationLim(5,2)) '$\Delta$'],...
            fig, [0.7625728553138 0.213743583907854 0.126685796465193 0.0295238091832116]);
    end
    hold
    
    %dC/dt for wave
    SSCChangeWave = dataWC(plotNum).waveDiffusionFlux+...
        dataWC(plotNum).CsedWave*dataWC(plotNum).ws+...
        dataWC(plotNum).vCsedWave;
    
    %dC/dt for current incase there is a miss balance
    SSCChangeSteady = dataWC(plotNum).Fdc+...
        dataWC(plotNum).CsedSteady*dataWC(plotNum).ws+...
        dataWC(plotNum).vCsedSteady;
    
    %Concentration change
    [AX,H1,H2] = plotIntegratedFluxFirst(SSCChangeWave/fluxNorm,...
        SSCChangeSteady/fluxNorm,...
        dataWC(plotNum).CsedSteady+dataWC(plotNum).CsedWave,y/delta,...
        bottomY,topY,plotSpecs,4,comp);
    
    %Plot diffusive Flux
    plotIntegratedFlux(dataWC(plotNum).waveDiffusionFlux/fluxNorm,...
        dataWC(plotNum).Fdc/fluxNorm,y/delta,bottomY,topY,...
        plotSpecs,1,comp,AX,H1,H2)
    
    %Plot settling Flux
    plotIntegratedFlux(dataWC(plotNum).CsedWave*dataWC(plotNum).ws/fluxNorm,...
        dataWC(plotNum).CsedSteady*dataWC(plotNum).ws/fluxNorm,y/delta,bottomY,topY,...
        plotSpecs,2,comp,AX,H1,H2)
    
    %Plot turb Flux
    plotIntegratedFlux(dataWC(plotNum).vCsedWave/fluxNorm,...
        dataWC(plotNum).vCsedSteady/fluxNorm,y/delta,bottomY,topY,...
        plotSpecs,3,comp,AX,H1,H2)
end
end



%% Functions
%%
function plotIntegratedFlux(Flux1,Flux2,z,z1,z2,plotSpecs,plotColorInd,comp,AX,H1,H2)
% This function computes the integrated sediment flux between z = z1 and z= z2
%
% Inputs:
% Flux1 - flux from currents
% Flux2 - flux from waves
% z - vector containing cell centered vertical coordinate
% z1,z2 - start and stop height of integration
% plotSpecs - specifications for plotting
% plotColorInd - index for the line color

indz1 = find(z>=z1,1,'first'); % find start index for integration
indz2 = find(z>=z2,1,'first'); % find end index for integration

[m,n] = size(Flux1);
steadyFlux = Flux2(indz2)-Flux2(indz1); % flux from currents into layer
if comp == 0
    plotVar = Flux1(indz2,:)-Flux1(indz1,:)+steadyFlux; % flux from waves into layer
elseif comp==1
    plotVar = Flux1(indz2,:)-Flux1(indz1,:); % flux from waves and currents
end

x = (0:n)/n*2*pi; %phase vector
if plotColorInd == 4
    plot(AX(1),x,[plotVar, plotVar(1)],...
        '-k','linewidth',plotSpecs.linewidth)
else
    plot(AX(1),x,[plotVar, plotVar(1)],...
        'color',plotSpecs.Cbudget(plotColorInd,:),...
        'linewidth',plotSpecs.linewidth)
end
end
%%
function [AX,H1,H2] = plotIntegratedFluxFirst(Flux1,Flux2,Csed,z,z1,z2,plotSpecs,plotColorInd,comp)
indz1 = find(z>=z1,1,'first');
indz2 = find(z>=z2,1,'first');

[m,n] = size(Flux1);
steadyFlux = Flux2(indz2)-Flux2(indz1);

if comp == 0
    plotVar = Flux1(indz2,:)-Flux1(indz1,:)+steadyFlux;
elseif comp==1
    plotVar = Flux1(indz2,:)-Flux1(indz1,:);
end

x = (0:n)/n*2*pi;
[AX,H1,H2]= plotyy(x,[plotVar, plotVar(1)],x,[Csed(1,:), Csed(1,1)]);

set(H1,'color','k','linewidth',plotSpecs.linewidth);
set(H2,'color','k','linewidth',plotSpecs.linewidth,'LineStyle','--');
set(AX,{'ycolor'},{'k';'k'})
set(AX(1),'ylim',[0 1],'FontSize',plotSpecs.axesfontsize)
set(AX(2),'ylim',[0 1],'FontSize',plotSpecs.axesfontsize)
set(gca,'xlim',xlim,'FontSize',plotSpecs.axesfontsize)
end

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
    'HorizontalAlignment','center',...
    'BackgroundColor',[1 1 1]);
end
