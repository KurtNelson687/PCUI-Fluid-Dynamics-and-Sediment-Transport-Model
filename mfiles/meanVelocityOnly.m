
%
% Filename    : post_process_pcui.m
% Author      : Goncalo Gil
% Description : Reads in PCUI data and plots a slice of the velocity field
%               and density field for the internal seiche test run.
%
% Author      : Goncalo Gil, Stanford University
% email       : gilg@stanford.edu
%

clear all; close all; clc;

isall = 0;
%working_folder = '/p/work1/knelson3/Waves4_1';
working_folder = '/p/work1/knelson3/ChannelFlow1Re2_7noEddy';
fname_xyz = 'xyz';
fname_uvw = 'output_UVW';

iskip = 1;
istart =  1;
iend = 12;%floor(params.nsteps/params.nsave)


% read the file containing the grid size and processor definitions
ftext = fileread(fullfile(working_folder, 'size.inc'));
params.ni = variable_value_pcui('ni',ftext);
params.nj = variable_value_pcui('nj',ftext);
params.nk = variable_value_pcui('nk',ftext);
params.px = variable_value_pcui('px',ftext);
params.py = variable_value_pcui('py',ftext);
params.pz = variable_value_pcui('pz',ftext);

% read the files containing the grid definition and assemble into array
% (does not include ghost grid points)
[~,y,~] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,isall,2);

'i have x,y,z data'


[m,n,p] = size(y);

for i=1:n
   if i == 1
      dy(i)= y(1,1,1)*2;
   else
      dy(i) = (y(1,i,1)-y(1,i-1,1)-dy(i-1)/2 )*2;
   end
end

count = 1;

u_depth_ave = zeros(length(istart:iskip:iend),1);
for istep = istart:iskip:iend
    istep 
    
    % Read velocity field
    [u,~,~] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
        params, 1,isall,2);

    u_mean_profile(count,:) = mean(mean(u,1),3);
   
    for i = 1:n
         u_depth_ave(count) = u_depth_ave(count)+u_mean_profile(count,i)*dy(i);
    end
    
    count = count+1;
end

u_depth_ave = u_depth_ave/(y(1,end,1)+dy(end)/2);

save([working_folder '/depthAveU.mat'],'u_mean_profile','u_depth_ave','dy')
