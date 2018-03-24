function [bedInd,viscousInd,bufferInd] = fixedZFluxes(dataWC,y,ydiff,plotSpecs,params,plotNum )
% Function: fixedZFluxes
%
% Author: Kurt Nelson
%
% Purpose: This function creates a 3-panel plot of 1) bed shear, 2) fluxes
% just above the bed, 3) fluxes at the top of the viscous sublayer layer,
% and 3) fluxes at the top of buffer layer.
%
% Inputs:
% 1) dataWC - PCUI data array created by data extractor that contains
%    velocity and sediment information
% 2) y - vector containing cell y coordinate
% 3) ydiff - vector containing vertical grid spacing
% 4) plotSpecs - array containing plot specifications
% 6) plotNum - run to plot
%
% Outputs:
% 1) bedInd - index of bottom most location for sediment fluxes
% 2) viscousInd - index of top of viscous sublayer
% 3) bufferInd - index of top of buffer layer
%%

%find index of top of buffer layer
ustar = sqrt(dataWC(plotNum).dpdxSteady*sum(ydiff)/params.rhoWater); % friction velocity
y_plus = y/(params.molecular_viscosity/ustar); % cell center in wall units
bedInd = 1; % bed index
viscousInd = find(y_plus>5,1,'first'); % index of top of viscous sublayer and also top of wave boundary layer
bufferInd = find(y_plus>30,1,'first'); % index of top of buffer layer

Cmean = 1/sum(ydiff)*ydiff*dataWC(plotNum).CsedSteady; % depth averaged sediment concentration for normalizing fluxes

fig = figure;
fig.PaperUnits = 'centimeters';
fig.PaperPosition =  [0 0 14 10];
clf
ha = tight_subplot(3,1,[0.03 .05],[.1 .125],[.1 .1]); %axes for plots

% x limits and labels for all plots
xlim = [0 2*pi];
ylim2 = [0 600];
ytick2 = [0, 150, 300, 450, 600];
xticklabels = {'0','\pi/4', '\pi/2','3\pi/4', '\pi','5\pi/4', '3\pi/2','7\pi/4', '2 \pi'};
xtick = [0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4,2*pi];

%% Panel 1 - bottom stress
axes('Position',[0.1 0.897619047619048 0.8 0.0895238095238177]) % position of top plot
set(gca,'box','on')
hold
temp = params.rhoWater*params.molecular_viscosity*dataWC(plotNum).umean(1,:)./y(1); % bed shear stress at each time step
phaseAveragedShear =postProcess(temp,0,y,0,'phaseAverage'); % phase average the shear stress
shearPhase = (0:length(phaseAveragedShear))/length(phaseAveragedShear)*2*pi; % phase of the shear stress

% compute erosion for flux at bed
index = abs(phaseAveragedShear)>=dataWC(plotNum).tauCrit;
dataWC(plotNum).E = zeros(size(phaseAveragedShear));
dataWC(plotNum).E(index) = dataWC(plotNum).Ased*(abs(phaseAveragedShear(index))-dataWC(plotNum).tauCrit);

plot(shearPhase,[phaseAveragedShear, phaseAveragedShear(1)]/dataWC(plotNum).tauCrit,...
    'linewidth',plotSpecs.linewidth)

plot(xlim,[-1 -1],':k','linewidth',1); % limit showing negative tau_crit
plot(xlim,[1 1],':k','linewidth',1); % limit showing positive tau_crit

myYlabel('$\frac{\langle \tau |_{z=0} \rangle_{\mathrm{p}}} {\tau_{crit}}$',plotSpecs.labfontsize)
set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel','','xlim',xlim)

%% Panel 2 - fluxes at the bed
axes(ha(1))
set(gca,'box','on')
set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel','')
hold
ylim1 = [-40 40];
ytick1 = [-40, -20, 0, 20, 40];
[~,~,~] = plotFixedZFluxes(dataWC,bedInd,plotSpecs,plotNum,Cmean,xlim,ylim1,ytick1,ylim2,ytick2); % plot sediment fluxes
addText('(a)', fig,[0.102831319877079 0.835238095408395 0.0328373602458409 0.0295238091832115]); % text for panel label

% y-label

% legend
leg = legend('$F_\mathrm{D}$','$F_\mathrm{s}$',...
    '$F_\mathrm{T}$','$\langle C \rangle_{\mathrm{p}}$');
set(leg,'interpreter','Latex','FontSize',plotSpecs.legfontsize,...
    'Position',[0.412756631328429 0.827572150242455 0.357672933169774 0.0307142853736877],...
    'Orientation','horizontal')

%% Panel 3 - fluxes at the top of the viscous sublayer
axes(ha(2))
set(gca,'box','on')
set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel','')
hold
ylim1 = [-12 12];
ytick1 = [-12, -6, 0, 6, 12];
[AX,~,~] = plotFixedZFluxes(dataWC,viscousInd,plotSpecs,plotNum,Cmean,xlim,ylim1,ytick1,ylim2,ytick2); % plot sediment fluxes
% y-label
ylabel(AX(1),'Normalized Sediment Fluxes','interpreter','Latex','FontSize',plotSpecs.labfontsize) % left y-axis
ylabel(AX(2),'$\langle C \rangle_{\mathrm{p}}$','interpreter','Latex','FontSize',plotSpecs.labfontsize) % left y-axis
addText('(b)', fig,[0.102136863606317 0.566190476360782 0.0342262727873666 0.0295238091832115]); % text for panel label

%% Panel 4 - fluxes at the top of the buffer layer
axes(ha(3))
set(gca,'box','on')
ylim1 = [-2.5 2.5];
ytick1 = [-2.5, -1.25, 0, 1.25, 2.5];
set(gca,'FontSize',plotSpecs.axesfontsize, 'xtick',xtick,'xticklabel',xticklabels)
hold
[~,~,~] = plotFixedZFluxes(dataWC,bufferInd,plotSpecs,plotNum,Cmean,xlim,ylim1,ytick1,ylim2,ytick2); % plot sediment fluxes
addText('(c)', fig,[0.0992598913056509 0.306666666836977 0.0328373602458409 0.0295238091832116]); % text for panel label

% x and y labels
myXlabel('\textbf{$ \theta $}',plotSpecs.labfontsize)
end

function [AX,H1,H2] = plotFixedZFluxes(dataWC,Ind,plotSpecs,plotNum,Cmean,xlim,ylim1,ytick1,ylim2,ytick2)
% This function plots unsteadyness, diffussion, settling, and turbluent
% sediment fluxes at a fixed height specified by
%
% Inputs:
% 1) dataWC - PCUI data array created by data extractor that contains
%    velocity and sediment information
% 2) y - vector containing cell y coordinate
% 3) bottomY - y value for bottom of the layer
% 4) topY - y value for the top of the layer
% 5) plotSpecs - array containing plot specifications
% 6) plotNum - run to plot

fluxNorm = dataWC(plotNum).ws*Cmean; % flux normalization
[~,n] = size(dataWC(plotNum).waveDiffusionFlux);
phase = (0:n)/n*2*pi; % phases to plot over

% Diffusive Flux
% combined wave and current diffusion flux
if Ind == 1 %if at the bed, specify the
    Combine_Fd = dataWC(plotNum).E;
else
    Combine_Fd = -1*(dataWC(plotNum).waveDiffusionFlux(Ind,:)+dataWC(plotNum).Fdc(Ind));
end
Combined_Csed = dataWC(plotNum).CsedWave(Ind,:)+dataWC(plotNum).CsedSteady(Ind);
[AX,H1,H2]= plotyy(phase,[Combine_Fd, Combine_Fd(1)]/fluxNorm,...
    phase,1000*[Combined_Csed, Combined_Csed(1)]);

% axes properties for plotyy
set(H1,'color',plotSpecs.Cbudget(1,:),'linewidth',plotSpecs.linewidth);
set(H2,'color','k','linewidth',plotSpecs.linewidth,'LineStyle','--');
set(AX,{'ycolor'},{'k';'k'})
set(AX(1),'xlim',xlim,'ylim',ylim1,'ytick',ytick1,'FontSize',plotSpecs.axesfontsize)
set(AX(2),'xlim',xlim,'ylim',ylim2,'ytick',ytick2,'FontSize',plotSpecs.axesfontsize)

% Settling Flux
% combined wave and current settling flux
Combine_Fs = -1*Combined_Csed *dataWC(plotNum).ws;
plot(AX(1),phase,[Combine_Fs, Combine_Fs(1)]/fluxNorm,'color',...
    plotSpecs.Cbudget(2,:),'linewidth',plotSpecs.linewidth)


% Turb Flux
% combined wave and current turbulent flux
Combine_Ft = -1*(dataWC(plotNum).vCsedWave(Ind,:)+dataWC(plotNum).vCsedSteady(Ind));
plot(AX(1),phase,[Combine_Ft, Combine_Ft(1)]/fluxNorm,'color',...
    plotSpecs.Cbudget(3,:),'linewidth',plotSpecs.linewidth)
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
    'HorizontalAlignment','center')
end


