

%
% Filename    : post_process_pcui.m
% Author      : Kurt Nelson
% Description : Reads in PCUI velocity data and computes dissipation
% email       : knelson@stanford.edu
%

clear all; close all; clc;

rhowater = 1000;
rhoSed = 2650;
istart =  1;
iend = 120;
iskip = 1;
isall = 0;
%working_folder = '/p/work1/knelson3/Waves4_1';
working_folder = '/work/knelson3/waveWithCurrents/dataWC200_11';
%working_folder = '/p/work1/knelson3/domainTest/Moser_3';
fname_xyz = 'xyz';
fname_Csed = 'output_Csed';
tic
% read the file containing the grid size and processor definitions
ftext = fileread(fullfile(working_folder, 'size.inc'));
params.ni = variable_value_pcui('ni',ftext);
params.nj = variable_value_pcui('nj',ftext);
params.nk = variable_value_pcui('nk',ftext);
params.px = variable_value_pcui('px',ftext);
params.py = variable_value_pcui('py',ftext);
params.pz = variable_value_pcui('pz',ftext);

count = 1;
for istep = istart:iskip:iend
   istep 
    Csed = read_binary_file_pcui(working_folder, fname_Csed, istep, ...
        params, 0,isall,2);  
    rho = rhowater+(1-rhowater/rhoSed).*Csed;
    [m,~,p] = size(rho);
    rhoProfile(:,count)=mean(mean(rho,1),3);
    rho = (rho-repmat(rhoProfile(:,count)',[m,1,p])).^2;
    rhoPrimeProfile(:,count)=mean(mean(rho,1),3);
    rhoPrimeProfile(:,count)=(rhoPrimeProfile(:,count)).^0.5;
count = count+1;
end
toc
save([working_folder '/meanRhoPrime.mat'],'rhoPrimeProfile','rhoProfile')
