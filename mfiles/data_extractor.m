
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

istart =  83;
isall = 1;
%working_folder = '/p/work1/knelson3/Waves4_1';
working_folder = '/p/work1/knelson3/ChannelFlow1_12';
fname_xyz = 'xyz';
fname_uvw = 'output_UVW';


% read the file containing the grid size and processor definitions
ftext = fileread(fullfile(working_folder, 'size.inc'));
params.ni = variable_value_pcui('ni',ftext);
params.nj = variable_value_pcui('nj',ftext);
params.nk = variable_value_pcui('nk',ftext);
params.px = variable_value_pcui('px',ftext);
params.py = variable_value_pcui('py',ftext);
params.pz = variable_value_pcui('pz',ftext);


for istep = istart
    istep 
    
    % Read velocity field
    [u,v,w] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
        params, 1,isall,2);

end

save([working_folder '/dataStep.mat'],'u','v','w')
