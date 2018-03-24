clear; close all;
data_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/QuadData/';
load('../DataFolder/CombinedHistData')
simulationType = 'histstrat200_';
dirnames = 1:15;
numRun = 1;

%Initialize
count = 0;
phase = 1;
for ndirs = 1:length(dirnames)
    if mod(numRun,2) == 0
        load([data_folder simulationType num2str(dirnames(ndirs)) 'cup'])
    else
        load([data_folder simulationType num2str(dirnames(ndirs))])
    end
    allHeightData(numRun).pdfHeights = pdfHeights;
    
    [numSteps,numHeights] = size(uvHisData);
    for i = 1:numHeights
        allEdgeData(numRun,i).uEdge = edgesAll(i).uEdge;
        allEdgeData(numRun,i).vEdge = edgesAll(i).vEdge;
        allEdgeData(numRun,i).rhoEdge = edgesAll(i).rhoEdge;
        for k = 1:numSteps
            if k < 13 && ndirs == 1
                allHistData(numRun,i,phase).uvMatrix = uvHisData(k,i).uvMatrix;
                allHistData(numRun,i,phase).vCsedMatrix = uvHisData(k,i).vCsedMatrix;
                %      allHistData(numRun,i,phase).vCsedMatrix = uvHisData(k,i).vCsedMatrix;
            else
                allHistData(numRun,i,phase).uvMatrix = allHistData(numRun,i,phase).uvMatrix + ...
                    uvHisData(k,i).uvMatrix;
                
                allHistData(numRun,i,phase).vCsedMatrix = allHistData(numRun,i,phase).vCsedMatrix + ...
                    uvHisData(k,i).vCsedMatrix;
                %  allHistData(numRun,i,phase).vCsedMatrix = allHistData(numRun,i,phase).vCsedMatrix + ...
                %      uvHisData(k,i).vCsedMatrix;
            end
            
            phase = phase+1;
            if phase==13
                phase = 1;
            end
        end
    end
    count = count+1;
end

save('../DataFolder/CombinedHistData','allHistData', 'allHeightData','allEdgeData')