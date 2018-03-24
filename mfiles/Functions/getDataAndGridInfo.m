% Script: getDataAndGridInfo
%
% Author: Kurt Nelson
%
% Purpose: Loads processed DNS data from PCUI, plot specification variables,
% grid data, and mean velocity profiles and turbulent intensities from 
% Moser et al. (1999) and del Alamo (2004).
% 
% Input: folder and file name for PCUI data
%%
% Load Del Alamo and Moser data 
load('/Users/kurtnelson/Desktop/Dropbox/Data/MoserData/MoserRe587') 
load('/Users/kurtnelson/Desktop/Dropbox/Data/DelAlamo/Re180')
load('/Users/kurtnelson/Desktop/Dropbox/Data/DelAlamo/Re180Profiles')
load('/Users/kurtnelson/Desktop/Dropbox/Data/DelAlamo/Re550')
load('/Users/kurtnelson/Desktop/Dropbox/Data/DelAlamo/Re550Profiles')

load([data_folder '/' data_fileWC]) % load data for waves and currents
load([data_folder '/' data_fileC]) % load data for currents only 
load([data_folder '/dataStepx']) % x coordinates from PCUI
load([data_folder '/dataStepy']) % y coordinates from PCUI
load([data_folder '/dataStepz']) % z coordinates from PCUI

load([data_folder '/figureProperties.mat']) %figure properties for plotting
load([data_folder '/paramsAndgridVal.mat']) %grid parameters

x = x(:,1,1); % make x a vector contianing only unique entries
dx = (x(2)-x(1)); % x grid spacing
z = z(1,1,:); % make z a vector contianing only unique entries
dz = z(2)-z(1); % z grid spacing
y = y(1,:,1); % make y a vector contianing only unique entries

% y grid spacing
dy(1) = 2*y(1);
for j = 2:length(y)
   dy(j) =  (y(j)-y(j-1)-dy(j-1)/2)*2;
end

H = y(end)+0.5*(y(end)-y(end-1)); % Height of domain 
