function [h] = TwoLayerPanel(dataWC,y,ydiff,plotSpecs,params,plotNum )
% Function: TwoLayerPanel
%
% Author: Kurt Nelson
%
% Purpose: This function creates a 5-panel plot of 1) pressure and
% free-stream velocity, 2) bed shear, 3) fluxes integrated from the bed to
% the top of the buffer layer, 4) fluxes integrated from the top of the buffer
% layer to the free-surface, 5) sediment conentration in each layer.
%
% Inputs:
% 1) dataWC - PCUI data array created by data extractor that contains
%    velocity and sediment information
% 2) y - vector containing cell y coordinate
% 3) ydiff - vector containing vertical grid spacing
% 4) plotSpecs - array containing plot specifications
% 5) plotNum - run to plot
%
%%

%find index of top of buffer layer
ustar = sqrt(dataWC(plotNum).dpdxSteady*sum(ydiff)/params.rhoWater); % friction velocity
y_plus = y/(params.molecular_viscosity/ustar); % cell center in wall units
yind = find(y_plus>30,1,'first'); % index of top of buffer layer
y_a = y(yind); % y value at the top of the buffer layer

fig = figure;
fig.PaperUnits = 'centimeters';
fig.PaperPosition =  plotSpecs.fullPageFigSize;
clf
ha = tight_subplot(5,1,[0.03 .05],[.1 .04],[.1 .04]);

% x limits and labels for all plots
xlim = [0 2*pi];
xticklabels = {'0','\pi/4', '\pi/2','3\pi/4', '\pi','5\pi/4', '3\pi/2','7\pi/4', '2 \pi'};
xtick = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4,2*pi];

%% Panel 1 - free-strem velocity and pressure
axes(ha(1))
set(gca,'box','on')
phase = 0:0.01:2*pi; % phase vector to plot over
plot(phase,sin(phase),'-',phase,cos(phase),'--')
set(gca,'xtick',xtick,'xticklabel',[],'xlim',[0 1],'ylim',[-1.1 1.1])

% y-label
myYlabel('$\frac{u_\mathrm{\infty}}{u_\mathrm{b}}$ and $\frac{S_\mathrm{W}}{u_\mathrm{b} \omega}$',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel','','xlim',xlim)

% legend
leg = legend('$u_\mathrm{\infty}$','$S_\mathrm{W}$');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.578568084750856 0.918571428741728 0.175006687641144 0.0307142853736877],...
    'Orientation','horizontal')

%% Panel 2 - bottom shear stress
axes(ha(2))
set(gca,'box','on')
hold
temp = params.rhoWater*params.molecular_viscosity*dataWC(plotNum).umean(1,:)./y(1); % bed shear stress at each time step
phaseAveragedShear =postProcess(temp,0,y,0,'phaseAverage'); % phase average the shear stress
shearPhase = (0:length(phaseAveragedShear))/length(phaseAveragedShear)*2*pi; % phase of the shear stress

plot(shearPhase,[phaseAveragedShear, phaseAveragedShear(1)]/dataWC(plotNum).tauCrit,...
    'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)

plot(xlim,[-1 -1],':k','linewidth',1); % limit showing negative tau_crit
plot(xlim,[1 1],':k','linewidth',1); % limit showing positive tau_crit

myYlabel('$\frac{\langle \tau |_{z=0} \rangle_{\mathrm{p}}} {\tau_{crit}}$',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel','','xlim',xlim)

%% Panel 3 - fluxes integrated from the bed to the top of the buffer layer
axes(ha(3))
set(gca,'box','on')
hold

Cmean = 1/sum(ydiff)*ydiff*dataWC(plotNum).CsedSteady; % depth averaged sediment concentration for normalizing fluxes
fluxNorm = dataWC(plotNum).ws*Cmean/(y_a-0); % flux normalization

plotAllFluxes(dataWC,y,y(1),y_a,plotSpecs,plotNum,fluxNorm) % plot sediment fluxes

% y-label
myYlabel('Sediment Fluxes',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel','','xlim',xlim)

% legend
leg = legend('$\frac{\partial C}{\partial t}$','$\Delta F_\mathrm{D,w}$',...
    '$\Delta F_\mathrm{w,w}$','$\Delta F_\mathrm{T,w}$');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.413688169331262 0.586763915914135 0.395095571449825 0.0313783730779376],...
    'Orientation','horizontal')

%% Panel 4 - fluxes integrated from top of the buffer layer to H
axes(ha(4))
set(gca,'box','on')
hold

fluxNorm = dataWC(plotNum).ws*Cmean/(y_a-0); % flux normalization
plotAllFluxes(dataWC,y,y_a,y(end),plotSpecs,plotNum,fluxNorm) % plot sediment fluxes

% y-label
myYlabel('Sediment Fluxes',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel','','xlim',xlim)

%% Panel 5 - sediment concentrations in the top and bottom layer
axes(ha(5))
set(gca,'box','on')
hold

phaseCsed =postProcess(dataWC(plotNum).Csed,0,y,0,'phaseAverage'); % phased-averaged sediment concentration
Clayer1 = 1/(y_a)*phaseCsed(1:yind,:)'*ydiff(1:yind)'; % layer-integrated phased-averaged sediment concentration in layer 1
Clayer2 = 1/(y(end)-y_a)*phaseCsed(yind:end,:)'*ydiff(yind:end)'; % layer-integrated phased-averaged sediment concentration in layer 2

plot(shearPhase,[Clayer1; Clayer1(1)]*1000,'linewidth',plotSpecs.linewidth)
plot(shearPhase,[Clayer2; Clayer2(1)]*1000,'linewidth',plotSpecs.linewidth)

% legend
leg = legend('bottom','top');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.597361184553762 0.224758297428597 0.199057727200644 0.0307142853736877],...
    'Orientation','horizontal')

% x and y labels
myXlabel('\textbf{$ \theta $}',plotSpecs.labfontsize)
myYlabel('SSC (mg/L)',plotSpecs.labfontsize)

set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel',xticklabels,'xlim',xlim)
end

function plotAllFluxes(dataWC,y,bottomY,topY,plotSpecs,plotNum,fluxNorm)
% This function plots unsteadyness, diffussion, settling, and turbluent
% sediment fluxes between the layer defined by bottomY and topY.
%
% Inputs:
% 1) dataWC - PCUI data array created by data extractor that contains
%    velocity and sediment information
% 2) y - vector containing cell y coordinate
% 3) bottomY - y value for bottom of the layer
% 4) topY - y value for the top of the layer
% 5) plotSpecs - array containing plot specifications
% 6) plotNum - run to plot

%dC/dt for wave
SSCChangeWave = dataWC(plotNum).waveDiffusionFlux+...
    dataWC(plotNum).CsedWave*dataWC(plotNum).ws+...
    dataWC(plotNum).vCsedWave;

plotIntegratedFlux(SSCChangeWave/fluxNorm,...
    dataWC(plotNum).Fdc/fluxNorm,y,bottomY,topY,...
    plotSpecs,4)

%Diffusive Flux
plotIntegratedFlux(dataWC(plotNum).waveDiffusionFlux/fluxNorm,...
    dataWC(plotNum).Fdc/fluxNorm,y,bottomY,topY,...
    plotSpecs,4)

%Diffusive Flux
plotIntegratedFlux(dataWC(plotNum).waveDiffusionFlux/fluxNorm,...
    dataWC(plotNum).Fdc/fluxNorm,y,bottomY,topY,...
    plotSpecs,1)

%Settling Flux
plotIntegratedFlux(dataWC(plotNum).CsedWave*dataWC(plotNum).ws/fluxNorm,...
    dataWC(plotNum).CsedSteady*dataWC(plotNum).ws/fluxNorm,y,bottomY,topY,...
    plotSpecs,2)

%Turb Flux
plotIntegratedFlux(dataWC(plotNum).vCsedWave/fluxNorm,...
    dataWC(plotNum).vCsedSteady/fluxNorm,y,bottomY,topY,...
    plotSpecs,3)
end

function plotIntegratedFlux(Flux1,Flux2,z,z1,z2,plotSpecs,plotColorInd)
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

[~,n] = size(Flux1);
steadyFlux = Flux2(indz2)-Flux2(indz1); % flux from currents into layer

plotVar = Flux1(indz2,:)-Flux1(indz1,:)+steadyFlux; % flux from waves and currents into the layer

x = (0:n)/n*2*pi; % phases to plot over
if plotColorInd == 4
    plot(x,[plotVar, plotVar(1)],...
        '-k','linewidth',plotSpecs.linewidth)
else
    plot(x,[plotVar, plotVar(1)],...
        'color',plotSpecs.Cbudget(plotColorInd,:),...
        'linewidth',plotSpecs.linewidth)
end
end

