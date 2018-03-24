close all; clear all; clear all;
addpath('./Functions');
data_folder = '/Users/kurtnelson/Desktop/Writing/WaveAndCurrent_showStratification/mfiles/Data';
%data_folder = '/Users/kurtnelson/Desktop/Writing/WaveAndCurrent_sedFocus/mfiles/Data';
save_folder = './binaryProfiles';
profile_name = 's2trat200cup';
run = 6;
StartPeriod = 1;
num_steps_period = 60;
aveStart = (StartPeriod-1)*num_steps_period;

%loading data
load([data_folder '/dataStrat'])
%load([data_folder '/dataWCfineSpinup'])
load([data_folder '/dataStepy']);
y = y(1,:,1);
H = y(end)+0.5*(y(end)-y(end-1));

for plotNum = run
    count = 1;
    [m,n] = size(dataWC(plotNum).umean);
    count = 1;
    for phase = 1:num_steps_period
        aveIndex = aveStart+phase:num_steps_period:n;
        umeanPhaseAve(:,count) = mean(dataWC(plotNum).umean(:,aveIndex),2);
        uvReyPhaseAve(:,count) = mean(dataWC(plotNum).uvRey(:,aveIndex),2);
        vCsedPhaseAve(:,count) = mean(dataWC(plotNum).vCsed(:,aveIndex),2);
        count = count+1;
    end
end

 urms = sqrt(mean(dataWC(plotNum).uTurb.^2,2));
 vrms = sqrt(mean(dataWC(plotNum).vTurb.^2,2));
rhorms = sqrt(mean(dataWC(plotNum).rhoPrimeMean.^2,2));


%% Plotting profiles
% Plotting average and all for particular phase
phase = 40;
figure
hold
aveIndex = aveStart+phase:num_steps_period:n;
cmap = colormap(cool(length(aveIndex)));
for i = 1:length(aveIndex)
   plot(dataWC(plotNum).umean(:,aveIndex(i)),y,'color',cmap(i,:))
end
plot(umeanPhaseAve(:,phase),y,'xk')
        
figure
hold
for phase = 1%:num_steps_period
        plot(umeanPhaseAve(:,phase),y)
end
load('1stStep.mat')
plot(umeanModel,yModel)
velDiff = umeanModel-umeanPhaseAve(:,phase);

save(['/Users/kurtnelson/Desktop/GeneralPCUIProcessingFiles/DataFolder/' profile_name],...
   'uvReyPhaseAve','vCsedPhaseAve')
% umeanPhaseAve = [-1*umeanPhaseAve(1,:); umeanPhaseAve; umeanPhaseAve(end,:)];
% write_binary_file_pcui_profile(save_folder, profile_name, params, umeanPhaseAve);
