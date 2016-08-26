%
% Filename    : initialize_pcui.m
% Author      : Goncalo Gil
% Description : Creates an input file to be read by PCUI. It can write
%               three- or four-dimensional arrays for the density and
%               velocity fields, respectively.
%
% Author      : Goncalo Gil, Stanford University
% email       : gilg@stanford.edu
%

clear all; close all; clc;
%working_folder = '/p/work1/knelson3/Waves1_1';
working_folder = '/work/knelson3/Domain4/T40_R1000_1';

sedConcentration = 0; %sediment concentration in mg/L
Ub = 0.1; %bottom orbital velocity

% These are the output files with the filenames stripped out of extensions
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
fname_rho = 'input_S';
fname_uvw = 'input_UVW';
fname_UVW_to_PCUI = 'uvw_init_from_matlab';
fname_rho_init_to_PCUI = 'rho_init_from_matlab';
fname_Csed_init_to_PCUI = 'Csed_init_from_matlab';
fname_rho_full_to_PCUI = 'rho_full_from_matlab';
fname_grid_to_PCUI = 'xyz_init_from_matlab';

% -------------------------------------------------------------------------
% Get problem parameters and variables from PCUI
% -------------------------------------------------------------------------
% read the file containing the parameter definition
ftext = fileread(fullfile(working_folder, 'pcuiRunParams.txt'));
params.dt = variable_value_pcui('dtime',ftext);
params.molecular_viscosity = variable_value_pcui('vis',ftext);
params.eddy_viscosity = variable_value_pcui('ak',ftext);
params.nsteps = variable_value_pcui('nstep',ftext);
params.nsave = variable_value_pcui('nsave',ftext);
params.rho_sed = variable_value_pcui('rhoSed',ftext);
params.rho_knot = variable_value_pcui('rhoWater',ftext);
params.dpdxSteady = variable_value_pcui('dpdxSteady',ftext);



% read the file containing the domain definition
ftext = fileread(fullfile(working_folder, 'cavity.f'));
params.wallUnit = variable_value_pcui('wallUnit',ftext);
params.numWallUnits = variable_value_pcui('numWallUnits',ftext);
params.maxRatio = variable_value_pcui('maxRatio',ftext);

% read the file containing the grid size and processor definitions
ftext = fileread(fullfile(working_folder, 'size.inc'));
params.ni = variable_value_pcui('ni',ftext);
params.nj = variable_value_pcui('nj',ftext);
params.nk = variable_value_pcui('nk',ftext);
params.px = variable_value_pcui('px',ftext);
params.py = variable_value_pcui('py',ftext);
params.pz = variable_value_pcui('pz',ftext);

params

%params.bx = params.wallUnit*params.numWallUnits*params.maxRatio*params.ni;
%params.by = variable_value_pcui('by',ftext);
%params.bz = variable_value_pcui('bz',ftext);

% read the files containing the grid definition and assemble into array
% (includes ghost grid points of each submatrix)
[x,y,z] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,1,2);

% -------------------------------------------------------------------------
% Initialize PCUI
% -------------------------------------------------------------------------
% Prepare density field
rho_init_pcui = ones(size(x))*1000;


% Initialize sediment concentration
%sedRadius = params.bx/5;
%sedx = params.bx/2;
%sedy = 4*params.by/5;
%sedz = params.bz/2;

[m,n,p] = size(x);
for i=1:m
    for j=1:n
        for k=1:p
            %                         Csed_init_pcui(i,j,k) = sedConcentration*10^-3;
            %                         rho_pert_pcui(i,j,k) = (1-params.rho_knot/params.rho_sed)*Csed_init_pcui(i,j,k);
           % temp = ((x(i,j,k)-sedx)^2+(y(i,j,k)-sedy)^2+(z(i,j,k)-sedz)^2)^0.5;
            %if temp <= sedRadius
                Csed_init_pcui(i,j,k) = sedConcentration*10^-3;
            %else
            %    Csed_init_pcui(i,j,k) = 0;
             %   rho_pert_pcui(i,j,k) = 0;
           %end
            
%                         if x(i,j,k)>=sedx-sedRadius && x(i,j,k)<=sedx + sedRadius && y(i,j,k)>=sedy-sedRadius && y(i,j,k)<=sedy + sedRadius && z(i,j,k)>=sedz-sedRadius && z(i,j,k)<=sedz + sedRadius
%                             Csed_init_pcui(i,j,k) = sedConcentration*10^-3;
%                             rho_pert_pcui(i,j,k) = (1-params.rho_knot/params.rho_sed)*Csed_init_pcui(i,j,k);
%                         else
%                             Csed_init_pcui(i,j,k) = 0;
%                             rho_pert_pcui(i,j,k) = 0;
%                         end
        end
    end
end
rho_full_pcui = rho_init_pcui+(1-params.rho_knot/params.rho_sed)*Csed_init_pcui;


% Write PCUI binary files depending on the number of processors
write_binary_file_pcui(working_folder, fname_rho_full_to_PCUI, params, rho_full_pcui);
write_binary_file_pcui(working_folder, fname_rho_init_to_PCUI, params, rho_init_pcui);
write_binary_file_pcui(working_folder, fname_Csed_init_to_PCUI, params, Csed_init_pcui);

clear rho_full_pcui rho_init_pcui Csed_init_pcui x z

%..........................................................................
%Initialize PCUI velocity
%..........................................................................
lasty = y(1,end-2,1)
u_star = sqrt(params.dpdxSteady*(y(1,end-2,1)+(y(1,end-2,1)-y(1,end-3,1))/2)/params.rho_knot)
zo = params.molecular_viscosity/(9*u_star);
u_log = zeros(size(y));
for i=1:m
    for j=3:n
        for k=1:p
	    if y(i,j,k)<=11*params.molecular_viscosity/u_star
		u_log(i,j,k) = y(i,j,k)/(params.molecular_viscosity/u_star)*u_star;
		else
		u_log(i,j,k) = u_star/0.41.*log(y(i,j,k)/zo);
	    end
        end
    end
end


max_u_log = Ub;%Use this line for wave init 
%max_u_log = max(max(max(u_log)))
a = -40/100*max_u_log;
b = 40/100*max_u_log;


uvw_pcui(:,:,:,1) = (b-a).*rand(size(y))+a;
%uvw_pcui(:,:,:,1) = u_log +(b-a).*rand(size(y))+a;
uvw_pcui(:,:,:,2) = (b-a).*rand(size(y))+a;
uvw_pcui(:,:,:,3) = (b-a).*rand(size(y))+a;

write_binary_file_pcui(working_folder, fname_UVW_to_PCUI, params, uvw_pcui);
