p = parpool('local',32);
spmd
     labindex
end
parfor a=1:32
   [a,labindex]
end
delete(p)
exit

