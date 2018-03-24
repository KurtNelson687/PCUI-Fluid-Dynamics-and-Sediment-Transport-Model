clear; close all;
data_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/QuadData/';
load('../DataFolder/CombinedQuadData')
simulationType = 'strat200_';
dirnames = 1:15;
numRun = 2;

%Initialize
count = 0;
phase = 1;
numphase = zeros(1,12);
for ndirs = 1:length(dirnames)
    if mod(numRun,2) == 0
        load([data_folder simulationType num2str(dirnames(ndirs)) 'cup'])
    else
        load([data_folder simulationType num2str(dirnames(ndirs))])
    end
    allHeightData(numRun).allHeights = allHeights;
    
    [numSteps,numHeights] = size(quadData);
    numThres = length(allThresholds);
    for i = 1:numHeights
        for j = 1:numThres
            for k = 1:numSteps
                if k < 13 && ndirs == 1
                    allQuadData(numRun,i,j).Nstress(:,phase) = quadData(k,i).Nstress(:,j);
                    allQuadData(numRun,i,j).Qstress(:,phase) = quadData(k,i).Qstress(:,j);
                    allQuadData(numRun,i,j).Nsed(:,phase) = quadData(k,i).Nsed(:,j);
                    allQuadData(numRun,i,j).Qsed(:,phase) = quadData(k,i).sed(:,j);
                    %if  j == 1
                     %   numphase(phase) = numphase(phase)+1;
                     %   vuratio(i,phase) = quadData(k,i).verticalAnisotropy;
                    %end
                else
                    allQuadData(numRun,i,j).Nstress(:,phase) = ...
                        allQuadData(numRun,i,j).Nstress(:,phase) + quadData(k,i).Nstress(:,j);
                    allQuadData(numRun,i,j).Qstress(:,phase) = ...
                        allQuadData(numRun,i,j).Qstress(:,phase) + quadData(k,i).Qstress(:,j);
                    allQuadData(numRun,i,j).Nsed(:,phase) = ...
                        allQuadData(numRun,i,j).Nsed(:,phase) + quadData(k,i).Nsed(:,j);
                    allQuadData(numRun,i,j).Qsed(:,phase) = ...
                        allQuadData(numRun,i,j).Qsed(:,phase) + quadData(k,i).sed(:,j);
                    %if  j == 1
                        %numphase(phase) = numphase(phase)+1;
                        %vuratio(i,phase) = quadData(k,i).verticalAnisotropy+ vuratio(i,phase);
                    %end
                end
                
                phase = phase+1;
                if phase==13
                    phase = 1;
                end
            end
        end
    end
    count = count+1;
end

save('../DataFolder/CombinedQuadData','allQuadData', 'allHeightData','allThresholds')