function [profile,demeanVar] = removeMean(var,type)
%%%%This function removes the horizontal mean of a variable
if strcmp(type,'step') == 1 
    [profile, demeanVar] = stepDemean(var);
end

end

function [outprofile, outDemean]= stepDemean(var)
[m,~,p,q] = size(var);
for l = 1:q
    outprofile(:,l)=mean(mean(var(:,:,:,l),1),3);
    outDemean(:,:,:,l) = var(:,:,:,l)-repmat(outprofile(:,l)',[m,1,p]);
end
end
