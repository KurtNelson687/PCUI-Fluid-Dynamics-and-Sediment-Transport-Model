function [count, sumVals] = quadAnalysis(x,y,xy,xyMean,allThresholds)
% Ratio of local value to mean
xyratio = abs(xy./xyMean);

% Find quadrant indicies
indQ1= find(x>0 & y>0);
indQ2= find(x<0 & y>0);
indQ3= find(x<0 & y<0);
indQ4= find(x>0 & y<0);

thresCount = 1;

for thres = allThresholds
    
    %Count and sum for quadrant 1
    xyratio1 = xyratio(indQ1);
    xy1 = xy(indQ1);
    indQ1thres = find(xyratio1>thres);
    count(1,thresCount)=length(indQ1thres);
    sumVals(1,thresCount)=sum(xy1(indQ1thres));
    
    %Count and sum for quadrant 2
    xyratio2 = xyratio(indQ2);
    xy2 = xy(indQ2);
    indQ2thres = find(xyratio2>thres);
    count(2,thresCount)=length(indQ2thres);
    sumVals(2,thresCount)=sum(xy2(indQ2thres));
    
    %Count and sum for quadrant 3
    xyratio3 = xyratio(indQ3);
    xy3 = xy(indQ3);
    indQ3thres = find(xyratio3>thres);
    count(3,thresCount)=length(indQ3thres);
    sumVals(3,thresCount)=sum(xy3(indQ3thres));
    
    %Count and sum for quadrant 4
    xyratio4 = xyratio(indQ4);
    xy4 = xy(indQ4);
    indQ4thres = find(xyratio4>thres);
    count(4,thresCount)=length(indQ4thres);
    sumVals(4,thresCount)=sum(xy4(indQ4thres));
    
    thresCount = thresCount + 1;
end

end