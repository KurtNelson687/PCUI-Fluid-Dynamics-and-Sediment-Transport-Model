clear; close all;
isall = 0;
addpath('../Functions');
addpath('./histcn')
load('../DataFolder/y.mat')
load('../DataFolder/strat200.mat')
load('../DataFolder/strat200VelRms.mat')
data_folder = '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/';
simulationType = 'StepTest1';
dirnames = 1;%[1 5 50 100 500 1000 1500 2000 5000];
fname_uvw = 'output_UVW';
fname_Csed = 'output_Csed';
iskip = 1;
istart =  1;
iend = 1;%floor(params.nsteps/params.nsave);
allThresholds = 0:0.5:10;
rhoWater = 1000;
pdfHeights = [1:10,15,20,25,30,35,40,35,50,60,70];
H = y(end)+(y(end)-y(end-1))/2;

for ndirs = 1:length(dirnames)
    working_folder = [data_folder simulationType];
    
    
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
    for istep = istart:iskip:iend
        % Read velocity field
        [u,v,~] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
            params, 1,isall,2);
        [m,n,p] = size(u);
        
        % Read Sediment field
        rho = read_binary_file_pcui(working_folder, fname_Csed, istep, ...
            params, 0,isall,2);
        rho = rhoWater+(1-rhoWater/params.rhoSed)*rho;
        
        heightCount = 1;
        %%%%%%%%%%%%%%
        heightCount2 = 1;
        %%%%%%%%%%%%%
        for height = allHeights
            Hindex = find(y/(params.molecular_viscosity/ustar) >= height,1,'first');
            
            % extract values and desired height
            % for Reynolds stress
            utemp = reshape(u(:,Hindex,:),[m*p,1]);
            utemp = utemp-mean(utemp);%umeanPhaseAve(Hindex,phase);
            vtemp = reshape(v(:,Hindex,:),[m*p,1]);
            uvtemp = utemp.*vtemp;
            
            % for vertical transport
            rhotemp = reshape(rho(:,Hindex,:),[m*p,1]);
            rhotemp = rhotemp-mean(rhotemp);
            rhovtemp = rhotemp.*vtemp;
            
            
            [quadData(stepCount,heightCount).Nstress,quadData(stepCount,heightCount).Qstress] = ...
                quadAnalysis(utemp,vtemp,uvtemp,uvReyPhaseAve(Hindex),allThresholds);
            
            [quadData(stepCount,heightCount).Nsed,quadData(stepCount,heightCount).sed] = ...
                quadAnalysis(rhotemp,vtemp,rhovtemp,vCsedPhaseAve(Hindex),allThresholds);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            valCheck = find(pdfHeights == height);
            if ~isempty(valCheck)
                binVals = [utemp,vtemp];
                edgesAll(heightCount2).uEdge = linspace(-5*urms(Hindex),5*urms(Hindex),200);
                edgesAll(heightCount2).vEdge = linspace(-5*vrms(Hindex),5*vrms(Hindex),200);
                [uvHisData(stepCount,heightCount2).uvMatrix,~,~,~] =...
                    histcn(binVals,edgesAll(heightCount2).uEdge,edgesAll(heightCount2).vEdge);
                heightCount2 = heightCount2 + 1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            heightCount = heightCount + 1;
        end
        
        stepCount = stepCount+1;
        if phase == 12
            phase = 1;
        else
            phase = phase+1;
        end
    end
    save([ working_folder '/quadData'],'quadData','allThresholds','allHeights')
    save([ working_folder '/histData'],'uvHisData','edgesAll','pdfHeights')
end