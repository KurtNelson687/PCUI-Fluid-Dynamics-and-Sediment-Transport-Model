function [outvar] = postProcess(var,x,y,z,type)
num_steps_period = 12*5;
if strcmp(type,'phaseAverage')==1
    [outvar]=phaseAverage(var,num_steps_period);
elseif strcmp(type,'depthAverage')==1
    [outvar]=depthAverage(var,y);
end
end

function [outvar]=phaseAverage(var,num_steps_period)
[~,endAve] =size(var);
for i = 1:num_steps_period
    aveIndex = i:num_steps_period:endAve;
    outvar(:,i) = mean(var(:,aveIndex),2)';
end
end

function [outvar]=depthAverage(var,y)
dy(1) = y(1)*2;
for j = 2:length(y)
    dy(j) = (y(j)-y(j-1)-dy(j-1)/2)*2;
end

[~,n] = size(var);

for i = 1:n
    outvar(i) = dy*var(:,i);
end
end
