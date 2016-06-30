%
% Filename    : post_process_pcui.m
% Author      : Goncalo Gil
% Description : Reads in PCUI data and plots a slice of the velocity field
%               and density field for the internal seiche test run.
%
% Author      : Goncalo Gil, Stanford University
% email       : gilg@stanford.edu
%

clear; close all; clc;

isall = 0;


%%%%%% Testing Simulations
%working_folder = '/home/knelson/Desktop/SettlingVelocityTest/10mgL';

dirnames = [1 5 50 100 500 1000 1500 2000 5000];
for ndirs = 1:length(dirnames)
working_folder = ['/home/knelson/Desktop/SettlingVelocityTest/' num2str(dirnames(ndirs)) 'mgL'];

% These are the output files with the filenames stripped out of extensions
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
fname_rho = 'output_rho';
fname_uvw = 'output_UVW';
fname_Csed = 'output_Csed';
fname_pressure = 'output_p';

% -------------------------------------------------------------------------
% Get problem parameters and variables from PCUI
% -------------------------------------------------------------------------
% read the file containing the parameter definition
ftext = fileread(fullfile(working_folder, 'io.f'));
params.dt = variable_value_pcui('dtime',ftext);
params.molecular_viscosity = variable_value_pcui('vis',ftext);
params.eddy_viscosity = variable_value_pcui('ak',ftext);
params.nsteps = variable_value_pcui('nstep',ftext);
params.nsave = variable_value_pcui('nsave',ftext);
params.dpdxSteady = variable_value_pcui('dpdxSteady',ftext);
params.waveMag = variable_value_pcui('waveMag',ftext);
params.Twave = variable_value_pcui('Twave',ftext);
params.waves = variable_value_pcui('waves',ftext);
params.ws = variable_value_pcui('ws',ftext);

% read the file containing the domain definition
ftext = fileread(fullfile(working_folder, 'cavity.f'));
params.bx = variable_value_pcui('bx',ftext);
params.by = variable_value_pcui('by',ftext);
params.bz = variable_value_pcui('bz',ftext);

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
[x,y,z] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,isall,2);


iskip = 1;
istart =  1;
iend = 60;%floor(params.nsteps/params.nsave);
tcount = 1;
for istep = istart:iskip:iend
    if ndirs==1
       time(tcount)=(istep-1)*params.dt*params.nsave; 
    end
    
    % Read velocity field
    [u,v,w] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
        params, 1,isall,2);
    
    % Read density field
    rho = read_binary_file_pcui(working_folder, fname_rho, istep, ...
        params, 0,isall,2);
  
    % Read Sediment field
    Csed = read_binary_file_pcui(working_folder, fname_Csed, istep, ...
        params, 0,isall,2);
    
    Xsed(tcount) = sum(sum(sum(x.*Csed)))/sum(sum(sum(Csed)));
    Ysed(tcount) = sum(sum(sum(y.*Csed)))/sum(sum(sum(Csed)));
    Zsed(tcount) = sum(sum(sum(z.*Csed)))/sum(sum(sum(Csed)));
    
    if tcount>1
       ws(tcount-1,ndirs)=sqrt((Xsed(tcount)-Xsed(tcount-1))^2+(Ysed(tcount)-Ysed(tcount-1))^2+(Zsed(tcount)-Zsed(tcount-1))^2)/(iskip*params.nsave*params.dt)/params.ws;
    end
   
    tcount=tcount+1;
end
clear Xsed Ysed Zsed
end

fig1 = figure(1);
hold
for i =1:length(dirnames)
   plot(time(2:end),ws(:,i)) 
end
xlabel('time (s)')
ylabel('Ws_{eff}/Ws_{grain}')


fig2 = figure(2);
plot(dirnames,mean(ws),'-k*') 
xlabel('C_{init} (mg/L)')
ylabel('Mean Ws_{eff}/Ws_{grain}')


