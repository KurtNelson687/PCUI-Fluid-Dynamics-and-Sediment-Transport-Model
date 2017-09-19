

%
% Filename    : post_process_pcui.m
% Author      : Kurt Nelson
% Description : Reads in PCUI velocity data and computes dissipation
% email       : knelson@stanford.edu
%

clear all; close all; clc;
round = 1;
if round == 1
istart =  1;
iend = 20;
elseif round == 2
istart =  21;
iend = 40;
elseif round == 3
istart =  41;
iend = 60;
elseif round == 4
istart =  61;
iend = 80;
elseif round == 5
istart =  81;
iend = 100;
elseif round == 6
istart =  101;
iend = 120;
end

iskip = 1;
isall = 0;
%working_folder = '/p/work1/knelson3/Waves4_1';
working_folder = '/work/knelson3/waveWithCurrents/dataWC350_8';
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
%%%% for i direction fft
temp = div(x,u(:,:,:),'i','fft');
dissFft = temp.^2;
temp = div(x,v(:,:,:),'i','fft');
dissFft = dissFft + temp.^2;
temp = div(x,w(:,:,:),'i','fft');
dissFft = dissFft + temp.^2;

%%%% for j direction 5 point
temp = div(y,u(:,:,:),'j','5point');
dissFft = dissFft+temp.^2;
temp = div(y,v(:,:,:),'j','5point');
dissFft = dissFft + temp.^2;
temp = div(y,w(:,:,:),'j','5point');
dissFft = dissFft + temp.^2;

%%%% for k direction
temp = div(z,u(:,:,:),'k','fft');
dissFft = dissFft+temp.^2;
temp = div(z,v(:,:,:),'k','fft');
dissFft = dissFft + temp.^2;
temp = div(z,w(:,:,:),'k','fft');
dissFft = dissFft + temp.^2;

dissFft = vis*dissFft;
meanDissFft(count,:) = mean(mean(dissFft,1),3);
clear dissFft 

count = count+1;
end
toc
yProfile = y(1,:,1);

if round == 1
save([working_folder '/meanDissFft1.mat'],'meanDissFft')
save([working_folder '/yProfile.mat'],'yProfile')
elseif round == 2
save([working_folder '/meanDissFft2.mat'],'meanDissFft')
elseif round == 3
save([working_folder '/meanDissFft3.mat'],'meanDissFft')
elseif round == 4
save([working_folder '/meanDissFft4.mat'],'meanDissFft')
elseif round == 5
save([working_folder '/meanDissFft5.mat'],'meanDissFft')
elseif round == 6
save([working_folder '/meanDissFft6.mat'],'meanDissFft')
end
