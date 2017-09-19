

%
% Filename    : post_process_pcui.m
% Author      : Kurt Nelson
% Description : Reads in PCUI velocity data and computes dissipation
% email       : knelson@stanford.edu
%

clear all; close all; clc;
istart =  1;
iend = 60;
iskip = 1;
isall = 0;
%working_folder = '/p/work1/knelson3/Waves4_1';
working_folder = '/work/knelson3/waveWithCurrents/dataWC200_1';
%working_folder = '/p/work1/knelson3/domainTest/Moser_3';
fname_xyz = 'xyz';
fname_uvw = 'output_UVW';
tic
vis = 10^-6;
% read the file containing the grid size and processor definitions
ftext = fileread(fullfile(working_folder, 'size.inc'));
params.ni = variable_value_pcui('ni',ftext);
params.nj = variable_value_pcui('nj',ftext);
params.nk = variable_value_pcui('nk',ftext);
params.px = variable_value_pcui('px',ftext);
params.py = variable_value_pcui('py',ftext);
params.pz = variable_value_pcui('pz',ftext);

[x,y,z] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,0,2);
count = 1;
for istep = istart:iskip:iend
    istep 
    % Read velocity field
    [utemp,v,w] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
        params, 1,isall,2);
[umean,u]= removeMean(utemp,'step');
clear utemp

parfor myid = 1:5
'im in'
%%%% for i direction fft
if myid == 1
temp1 = div(x,u(:,:,:),'i','fft');
elseif myid == 2
temp2 = div(y,w(:,:,:),'j','5point');
elseif myid == 3
temp3 = div(y,v(:,:,:),'j','5point');
elseif myid == 4
temp4 = div(y,u(:,:,:),'j','5point');
elseif myid == 5
temp5 = div(z,u(:,:,:),'k','fft');
end
end
'through u and y'
clear u 
dissFft =  temp1.^2+temp2.^2+temp3.^2+temp4.^2+temp5.^2;
clear temp1 temp2 temp3 temp4 temp5


parfor myid = 1:4
%%%% for i direction fft
if myid == 1
temp1 = div(x,v(:,:,:),'i','fft');
elseif myid == 2
temp2 = div(x,w(:,:,:),'i','fft');
elseif myid == 3
temp3 = div(z,v(:,:,:),'k','fft');
elseif myid == 4
temp4 = div(z,w(:,:,:),'k','fft');
end
end

dissFft =  dissFft+temp1.^2+temp2.^2+temp3.^2+temp4.^2;
clear temp1 temp2 temp3 temp4

dissFft = vis*dissFft;
meanDissFft(count,:) = mean(mean(dissFft,1),3);
count = count+1;
end
toc
yProfile = y(1,:,1);
save([working_folder '/yProfile.mat'],'yProfile')
save([working_folder '/meanDissFft1.mat'],'meanDissFft')
