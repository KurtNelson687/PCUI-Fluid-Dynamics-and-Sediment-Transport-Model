clear all ; close all;

dirnames = 1:15;
addpath('./histcn')
poolobj = parpool('local',length(dirnames))
runNum = 1;
isall = 0;
load('~/dataForQuad/y.mat')
data_folder = '/p/work1/knelson3/waveWithCurrents/';
fname_uvw = 'output_UVW';
fname_Csed = 'output_Csed';
iskip = 1;
istart =  1;
iend = 120;%floor(params.nsteps/params.nsave);

H = y(end)+(y(end)-y(end-1))/2;

if runNum == 1
    simulationType = 'strat200_';
    load('~/dataForQuad/strat200.mat')
    load('~/dataForQuad/strat200cupVelRms.mat')
elseif runNum == 2
    simulationType = 'strat200_';
    load('~/dataForQuad/strat200cup.mat')
    load('~/dataForQuad/strat200cupVelRms.mat')
elseif runNum == 3
    simulationType = 'strat350_';
    load('~/dataForQuad/strat350.mat')
    load('~/dataForQuad/strat350cupVelRms.mat')
elseif runNum == 4
    simulationType = 'strat350_';
    load('~/dataForQuad/strat350cup.mat')
    load('~/dataForQuad/strat350cupVelRms.mat')
elseif runNum == 5
    simulationType = 'strat500_';
    load('~/dataForQuad/strat500.mat')
    load('~/dataForQuad/strat500cupVelRms.mat')
elseif runNum == 6
    simulationType = 'strat500_';
    load('~/dataForQuad/strat500cup.mat')
    load('~/dataForQuad/strat500cupVelRms.mat')
end

spmd(length(dirnames))
    if mod(runNum,2) == 0
        working_folder = [data_folder simulationType num2str(dirnames(labindex)) 'cup'];
    else
        working_folder = [data_folder simulationType num2str(dirnames(labindex))];
    end
    
    
    % -------------------------------------------------------------------------
    % Get problem parameters and variables from PCUI
    % -------------------------------------------------------------------------
    % read the file containing the parameter definition
    ftext = fileread(fullfile(working_folder, 'pcuiRunParams.txt'));
    params.molecular_viscosity = variable_value_pcui('vis',ftext);
    params.dpdxSteady = variable_value_pcui('dpdxSteady',ftext);
    params.rhoSed = variable_value_pcui('rhoSed',ftext);
    ustar = sqrt(params.dpdxSteady*H/rhoWater);
    allHeights = [1:60,70:20:max(y/(params.molecular_viscosity/ustar))];
    
    
    % read the file containing the grid size and processor definitions
    ftext = fileread(fullfile(working_folder, 'size.inc'));
    params.ni = variable_value_pcui('ni',ftext);
    params.nj = variable_value_pcui('nj',ftext);
    params.nk = variable_value_pcui('nk',ftext);
    params.px = variable_value_pcui('px',ftext);
    params.py = variable_value_pcui('py',ftext);
    params.pz = variable_value_pcui('pz',ftext);
    
    stepCount = 1;
    heightCount = 1;
    phase = 1;
    
    saveQuad(working_folder,quadData,allThresholds,allHeights,uvHisData,edgesAll,pdfHeights)
end

delete(poolobj)
