%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The functions found in this file create differnt plots
% from PCUI data. The main function is called
% by "HPC_movie.m".
%
%
% Author      : Kurt Nelson, Stanford University
% email       : knelson3@stanford.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WCplottingFunctions(x_plotvar1,y_plotvar1,plot_type,plotNum,params,plotSpecs,xplotNorm,yplotNorm,phase)
num_steps_period = 12*5;
if strcmp(plot_type,'timePlot')==1
    plotTimeSeries(plotNum,x_plotvar1,y_plotvar1,params,plotSpecs,xplotNorm,yplotNorm)
elseif strcmp(plot_type,'WaveTheo')==1
    plotTheoWave(phase,y_plotvar1,params,xplotNorm,yplotNorm)
elseif strcmp(plot_type,'periodAvePlot')==1
    plotPeriodAve(plotNum,x_plotvar1,y_plotvar1,params,plotSpecs,xplotNorm,yplotNorm,num_steps_period)
elseif strcmp(plot_type,'phaseAvePlot')==1
    plotPhaseAve(plotNum,x_plotvar1,y_plotvar1,params,plotSpecs,xplotNorm,yplotNorm,phase,num_steps_period)
elseif strcmp(plot_type,'phaseAvePlotTurb')==1
    plotPhaseAveTurb(plotNum,x_plotvar1,y_plotvar1,params,plotSpecs,xplotNorm,yplotNorm,phase,num_steps_period)
elseif strcmp(plot_type,'phaseAvePlotLogx')==1
    plotPhaseAveLogx(plotNum,x_plotvar1,y_plotvar1,params,plotSpecs,xplotNorm,yplotNorm,phase,num_steps_period)
elseif strcmp(plot_type,'timeAveProfilePlot')==1
    plotTimeAveProfile(plotNum,x_plotvar1,y_plotvar1,params,xplotNorm,yplotNorm,plotSpecs,phase,num_steps_period)
elseif strcmp(plot_type,'timeAveProfilePlotTurb')==1
    plotTimeAveProfileTurb(plotNum,x_plotvar1,y_plotvar1,params,xplotNorm,yplotNorm,plotSpecs,phase,num_steps_period)
elseif strcmp(plot_type,'LumleyDiagram')==1
    plotLumleyDiagram
elseif strcmp(plot_type,'Lumley')==1
    plotLumley(x_plotvar1,y_plotvar1,plotNum,plotSpecs)
 elseif strcmp(plot_type,'WaveComPhase')==1
    plotWaveComponentPhase(plotNum,x_plotvar1,y_plotvar1,xplotNorm,yplotNorm,plotSpecs,phase,num_steps_period)
end
end

function plotTimeSeries(plotNum,x_plotvar1,y_plotvar1,params,plotSpecs,xplotNorm,yplotNorm)
plot(x_plotvar1/xplotNorm,y_plotvar1/yplotNorm,plotSpecs.linetype{plotNum},'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
end

function plotPeriodAve(plotNum,x_plotvar1,y_plotvar1,params,plotSpecs,xplotNorm,yplotNorm, num_steps_period)
numAve = floor(length(x_plotvar1)/ num_steps_period);
for i = 1:numAve
    index = (i-1)* num_steps_period+1;
    xplot(i) = mean(x_plotvar1(index:index+ num_steps_period-1));
    yplot(i) = mean(y_plotvar1(index:index+ num_steps_period-1));
end
plot(xplot/xplotNorm,yplot/yplotNorm,plotSpecs.linetype{plotNum},'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
end

function plotPhaseAve(plotNum,x_plotvar1,y_plotvar1,params,plotSpecs,xplotNorm,yplotNorm,phase,num_steps_period)
[~,endAve] =size(x_plotvar1);
aveIndex = phase:num_steps_period:endAve;
xplot = mean(x_plotvar1(:,aveIndex),2)';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% endAve = floor(endAve/num_steps_period)*num_steps_period;
% aveIndex = 1:endAve;
% 
% %%This is compute the differnce between time avearge and phase average profiles
% endAvePeriod = floor(endAve/num_steps_period)*num_steps_period;
% aveIndexPeriod = 1:endAvePeriod;
% xplotPeriod = mean(x_plotvar1(:,aveIndexPeriod),2)';
% 
% xplot = abs(xplot-xplotPeriod)./xplot*100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%varXplot = mean(var(x_plotvar1(1:15,aveIndex)-xplot(1:15)'))
plot(xplot/xplotNorm,y_plotvar1/yplotNorm,plotSpecs.linetype{plotNum},'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
end

function plotPhaseAveTurb(plotNum,x_plotvar1,y_plotvar1,params,plotSpecs,xplotNorm,yplotNorm,phase,num_steps_period)
[~,endAve] =size(x_plotvar1);
aveIndex = phase:num_steps_period:endAve;
xplot = sqrt(mean(x_plotvar1(:,aveIndex).^2,2))';
%varXplot = mean(var(x_plotvar1(1:15,aveIndex)-xplot(1:15)'))
plot(xplot/xplotNorm,y_plotvar1/yplotNorm,plotSpecs.linetype{plotNum},'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
end

function plotPhaseAveLogx(plotNum,x_plotvar1,y_plotvar1,params,plotSpecs,xplotNorm,yplotNorm,phase,num_steps_period)
[~,endAve] =size(x_plotvar1);
aveIndex = phase:num_steps_period:endAve;
xplot = mean(x_plotvar1(:,aveIndex),2)';

semilogx(xplot/xplotNorm,y_plotvar1/yplotNorm,plotSpecs.linetype{plotNum},'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
end

function plotTheoWave(phase,y_plotvar1,params,xplotNorm,yplotNorm)
u_theo = cos(phase)-exp(-y_plotvar1/params.delta).*cos(phase-y_plotvar1/params.delta);
plot(u_theo/xplotNorm,y_plotvar1/yplotNorm,'xk','markersize', 4)
end

function plotTimeAveProfile(plotNum,x_plotvar1,y_plotvar1,params,xplotNorm,yplotNorm,plotSpecs,phase,num_steps_period)
[~,endAve] =size(x_plotvar1);
endAve = floor(endAve/num_steps_period)*num_steps_period;
aveIndex = 1:endAve;

% startPeriod = 60;
% aveIndex = (startPeriod-1)*60+1:endAve;

% startPeriod = 1;
% endPeriod = 90;
% aveIndex = (startPeriod-1)*60+1:endPeriod*60;

xplot = mean(x_plotvar1(:,aveIndex),2)';
ustar = sqrt(params.molecular_viscosity*xplot(1)/y_plotvar1(1))

if phase == 0
    plot(xplot/xplotNorm,y_plotvar1/yplotNorm,plotSpecs.linetype{plotNum},'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
elseif phase == 1
    semilogx(y_plotvar1/(params.molecular_viscosity/ustar),xplot/ustar,plotSpecs.linetype{plotNum},'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
end
end

function plotTimeAveProfileTurb(plotNum,x_plotvar1,y_plotvar1,params,xplotNorm,yplotNorm,plotSpecs,phase,num_steps_period)
[~,endAve] =size(x_plotvar1);
endAve = floor(endAve/num_steps_period)*num_steps_period;
aveIndex = 1:endAve;

%startPeriod = 60;
%aveIndex = (startPeriod-1)*60+1:endAve;

% startPeriod = 1;
% endPeriod = 90;
% aveIndex = (startPeriod-1)*60+1:endPeriod*60;

xplot = sqrt(mean(x_plotvar1(:,aveIndex).^2,2))';
ustar = sqrt(params.molecular_viscosity*xplot(1)/y_plotvar1(1));

if phase == 0
    plot(xplot/xplotNorm,y_plotvar1/yplotNorm,plotSpecs.linetype{plotNum},'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
elseif phase == 1
    semilogx(y_plotvar1/(params.molecular_viscosity/ustar),xplot/ustar,plotSpecs.linetype{plotNum},'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
end
end

function plotLumleyDiagram
x1 = -0.01:0.0001:2/27;
y1 = 27/9*x1+1/9;
plot(x1,y1,':k','linewidth',0.5)

x2 = 0:0.0001:2/27;
y2 = 3*(1/2*x2).^(2/3);
plot(x2,y2,':k','linewidth',0.5)

x3 = -0.01:0.0001:0;
y3 = 3*(-1/2*x3).^(2/3);
plot(x3,y3,':k','linewidth',0.5)
end

function plotLumley(x_plotvar1,y_plotvar1,plotNum,plotSpecs)
 plot(x_plotvar1,y_plotvar1,plotSpecs.linetype{plotNum},'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
end

function plotWaveComponentPhase(plotNum,x_plotvar1,y_plotvar1,xplotNorm,yplotNorm,plotSpecs,phase,num_steps_period)
[~,endAve] =size(x_plotvar1);
endAve = floor(endAve/num_steps_period)*num_steps_period;
aveIndex = 1:endAve;

meanProfile = mean(x_plotvar1(:,aveIndex),2)';

[~,endAve] =size(x_plotvar1);
aveIndex = phase:num_steps_period:endAve;
xplot = mean(x_plotvar1(:,aveIndex),2)';

xplot = xplot-meanProfile;
plot(xplot/xplotNorm,y_plotvar1/yplotNorm,plotSpecs.linetype{plotNum},'color',plotSpecs.C(plotNum,:),'linewidth',plotSpecs.linewidth)
end
