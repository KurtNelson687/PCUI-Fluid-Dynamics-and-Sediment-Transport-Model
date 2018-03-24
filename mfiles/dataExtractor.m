close all; clear all;
addpath('./Functions');
append = 1:50;
replace_profile = 6;
iskip = 1;
include_post = 0;
foldCount = 1;
for i = append
    %folderName{foldCount} = ['strat500_' num2str(i)];
    folderName{foldCount} = ['strat500_' num2str(i) 'cup'];
    dataFolder{foldCount}= '/Users/kurtnelson/Desktop/HPC_transfer/StratStudy/';
    foldCount = foldCount+1;
end

working_folder = [dataFolder{1} folderName{1}];
% These are the output files with the filenames stripped out of extensions
% (extensions are chosen automatically based on number of processors).
fname_time = 'outputp_time';
fname_time_full = 'output_time';
fname_kount= 'outputpVal_kount';
fname_umean = 'outputp_umean';
fname_vmean = 'outputp_vmean';
fname_wmean = 'outputp_wmean';
fname_uTurb = 'outputp_uTurb';
fname_vTurb = 'outputp_vTurb';
fname_wTurb = 'outputp_wTurb';
fname_uvRey = 'outputp_uvRey';
fname_uwRey = 'outputp_uwRey';
fname_vwRey = 'outputp_vwRey';
fname_DisDepth = 'outputpval_DisDepth';
fname_ProDepth  = 'outputpval_ProDepth';
fname_phase  = 'outputpval_wavephase';

fname_vstMean = 'outputp_vstMean';
fname_rruMean = 'outputp_rruMean';
fname_uDepthAve = 'outputpVal_udepth';
fname_Source = 'outputpVal_Source';
fname_drive = 'outputpVal_drive';
fname_kinetic = 'outputpVal_kineticdepth';
fname_dissipation = 'outputp_dissmean';
fname_cfl = 'outputpval_cfl';
fname_Sed = 'outputp_cSed';
fname_vCSed = 'outputp_vCSed';
fname_BruntN = 'outputp_BruntN';
fname_rhoMean = 'outputp_rhoMean';
fname_rhoSqrMean = 'outputp_rhoSqrMean';
fname_rhoPrimeMean = 'outputp_rhoPrimeMean';
fname_sedTotal1 = 'outputpVal_sedTotal1';
fname_sedTotal2 = 'outputpVal_sedTotal2';

% read the file containing the grid size and processor definitions
ftext = fileread(fullfile(working_folder, 'size.inc'));
params.ni = variable_value_pcui('ni',ftext);
params.nj = variable_value_pcui('nj',ftext);
params.nk = variable_value_pcui('nk',ftext);
params.px = variable_value_pcui('px',ftext);
params.py = variable_value_pcui('py',ftext);
params.pz = variable_value_pcui('pz',ftext);


ftext = fileread(fullfile(working_folder, 'pcuiRunParams.txt'));
Add_data(1).dt = variable_value_pcui('dtime',ftext);
Add_data(1).nsave = variable_value_pcui('nsave',ftext);
Add_data(1).nsavePro = variable_value_pcui('nsavePro',ftext);
Add_data(1).dpdxSteady = variable_value_pcui('dpdxSteady',ftext);
Add_data(1).dpdxWave = variable_value_pcui('dpdxWave',ftext);
Add_data(1).iTKE = variable_value_pcui('iTKE',ftext);
Add_data(1).ieddy = variable_value_pcui('ieddy',ftext);
Add_data(1).waves = variable_value_pcui('waves',ftext);
Add_data(1).Twave = variable_value_pcui('Twave',ftext);

Add_data(1).ised = variable_value_pcui('ised',ftext);
if Add_data(1).ised == 1
    Add_data(1).rhoSed = variable_value_pcui('rhoSed',ftext);
    Add_data(1).tauCrit = variable_value_pcui('tauCrit',ftext);
    Add_data(1).Ased = variable_value_pcui('Ased',ftext);
    Add_data(1).nsed = variable_value_pcui('nsed',ftext);
    Add_data(1).ws = variable_value_pcui('ws',ftext);
end
Add_data(1).ak = variable_value_pcui('ak',ftext);
params.rhoWater = variable_value_pcui('rhoWater',ftext);
params.molecular_viscosity = variable_value_pcui('vis',ftext);

array_index = 1;

f = dir([working_folder '/' fname_wmean '.100']);
istop = f.bytes/(8.1250*params.nj/params.py);
count = 1;

for istep = 1:iskip:istop
    Add_data(array_index).umean(:,count) = read_binary_file_pcui_profile(working_folder, fname_umean, istep,params,0,0);
    Add_data(array_index).vmean(:,count) = read_binary_file_pcui_profile(working_folder, fname_vmean, istep,params,0,0);
    Add_data(array_index).wmean(:,count) = read_binary_file_pcui_profile(working_folder, fname_wmean, istep,params,0,0);
    Add_data(array_index).uTurb(:,count) = read_binary_file_pcui_profile(working_folder, fname_uTurb, istep,params,0,0);
    Add_data(array_index).vTurb(:,count) = read_binary_file_pcui_profile(working_folder, fname_vTurb, istep,params,0,0);
    Add_data(array_index).wTurb(:,count) = read_binary_file_pcui_profile(working_folder, fname_wTurb, istep,params,0,0);
    Add_data(array_index).uvRey(:,count) = read_binary_file_pcui_profile(working_folder, fname_uvRey, istep,params,0,0);
    Add_data(array_index).uwRey(:,count) = read_binary_file_pcui_profile(working_folder, fname_uwRey, istep,params,0,0);
    Add_data(array_index).vwRey(:,count) = read_binary_file_pcui_profile(working_folder, fname_vwRey, istep,params,0,0);
    Add_data(array_index).dissipation(:,count) = read_binary_file_pcui_profile(working_folder, fname_dissipation, istep,params,0,0);
    Add_data(array_index).uDepthAve(1,count) = read_binary_file_pcui_value(working_folder, fname_uDepthAve, istep,1);
    Add_data(array_index).time(1,count) = read_binary_file_pcui_value(working_folder, fname_time, istep,1);
    Add_data(array_index).kount(1,count) = read_binary_file_pcui_value(working_folder, fname_kount, istep,1);
    Add_data(array_index).Sdrive(1,count) = read_binary_file_pcui_value(working_folder, fname_drive, istep,1);
    Add_data(array_index).totalKinetic(1,count) = read_binary_file_pcui_value(working_folder, fname_kinetic, istep,1);
    %Add_data(array_index).cfl(1,count) = read_binary_file_pcui_value(working_folder,   fname_cfl, istep,1);
    if Add_data(1).iTKE == 1
        Add_data(array_index).proDepth(1,count) = read_binary_file_pcui_value(working_folder,   fname_ProDepth, istep,1);
        Add_data(array_index).disDepth(1,count) = read_binary_file_pcui_value(working_folder,   fname_DisDepth, istep,1);
    end
    if Add_data(1).waves == 1
        Add_data(array_index).phase(1,count) = read_binary_file_pcui_value(working_folder,   fname_phase, istep,1);
    end
    
    if Add_data(1).ised == 1
        Add_data(array_index).Csed(:,count) = read_binary_file_pcui_profile(working_folder, fname_Sed, istep,params,0,0);
        Add_data(array_index).vCsed(:,count) = read_binary_file_pcui_profile(working_folder, fname_vCSed, istep,params,0,0);
        Add_data(array_index).CsedTotal1(1,count) = read_binary_file_pcui_value(working_folder,   fname_sedTotal1, istep,1);
        Add_data(array_index).CsedTotal2(1,count) = read_binary_file_pcui_value(working_folder,   fname_sedTotal2, istep,1);
        Add_data(array_index).BruntN(:,count) = read_binary_file_pcui_profile(working_folder, fname_BruntN, istep,params,0,0);
        Add_data(array_index).rhoMean(:,count) = read_binary_file_pcui_profile(working_folder, fname_rhoMean, istep,params,0,0);
        Add_data(array_index).rhoSqrMean(:,count) = read_binary_file_pcui_profile(working_folder, fname_rhoSqrMean, istep,params,0,0);
        Add_data(array_index).rhoPrimeMean(:,count) = read_binary_file_pcui_profile(working_folder, fname_rhoPrimeMean, istep,params,0,0);
    end
    
    if Add_data(1).ieddy == 1
        Add_data(array_index).vst(:,count) = read_binary_file_pcui_profile(working_folder, fname_vstMean, istep,params,0,0);
        Add_data(array_index).rr(:,count) = read_binary_file_pcui_profile(working_folder, fname_rruMean, istep,params,0,0);
    end
    count =count+1;
end

allFolder = folderName{1};


count
if append == 1
    Add_data(array_index).steps = istep;
else
    
    for numFolder = 2:length(append)
        allFolder  = [allFolder  ' and ' folderName{numFolder}];
        working_folder = [dataFolder{numFolder} folderName{numFolder}];
        f = dir([working_folder '/' fname_wmean '.100']);
        istop = f.bytes/(8.1250*params.nj/params.py);
        
        temp = read_binary_file_pcui_value(working_folder, fname_time, 1,1);
        count = find(Add_data(1).time==temp,1,'first');
        Add_data(array_index).steps = istop+count;
        if isempty(count)
            'no time match'
            count = length(Add_data.time)+1;
            Add_data(array_index).steps = istop+count-1;
        end
        
        
        for istep = 1:iskip:istop
            Add_data(array_index).umean(:,count) = read_binary_file_pcui_profile(working_folder, fname_umean, istep,params,0,0);
            Add_data(array_index).vmean(:,count) = read_binary_file_pcui_profile(working_folder, fname_vmean, istep,params,0,0);
            Add_data(array_index).wmean(:,count) = read_binary_file_pcui_profile(working_folder, fname_wmean, istep,params,0,0);
            Add_data(array_index).uTurb(:,count) = read_binary_file_pcui_profile(working_folder, fname_uTurb, istep,params,0,0);
            Add_data(array_index).vTurb(:,count) = read_binary_file_pcui_profile(working_folder, fname_vTurb, istep,params,0,0);
            Add_data(array_index).wTurb(:,count) = read_binary_file_pcui_profile(working_folder, fname_wTurb, istep,params,0,0);
            Add_data(array_index).uvRey(:,count) = read_binary_file_pcui_profile(working_folder, fname_uvRey, istep,params,0,0);
            Add_data(array_index).uwRey(:,count) = read_binary_file_pcui_profile(working_folder, fname_uwRey, istep,params,0,0);
            Add_data(array_index).vwRey(:,count) = read_binary_file_pcui_profile(working_folder, fname_vwRey, istep,params,0,0);
            Add_data(array_index).dissipation(:,count) = read_binary_file_pcui_profile(working_folder, fname_dissipation, istep,params,0,0);
            Add_data(array_index).time(1,count) = read_binary_file_pcui_value(working_folder, fname_time, istep,1);
            Add_data(array_index).kount(1,count) = read_binary_file_pcui_value(working_folder, fname_kount, istep,1);
            Add_data(array_index).uDepthAve(1,count) = read_binary_file_pcui_value(working_folder, fname_uDepthAve, istep,1);
            Add_data(array_index).Sdrive(1,count) = read_binary_file_pcui_value(working_folder, fname_drive, istep,1);
            Add_data(array_index).totalKinetic(1,count) = read_binary_file_pcui_value(working_folder, fname_kinetic, istep,1);
            %   Add_data(array_index).cfl(1,count) = read_binary_file_pcui_value(working_folder,   fname_cfl, istep,1);
            if Add_data(1).iTKE == 1
                Add_data(array_index).proDepth(1,count) = read_binary_file_pcui_value(working_folder,   fname_ProDepth, istep,1);
                Add_data(array_index).disDepth(1,count) = read_binary_file_pcui_value(working_folder,   fname_DisDepth, istep,1);
            end
            if Add_data(1).waves == 1
                Add_data(array_index).phase(1,count) = read_binary_file_pcui_value(working_folder,   fname_phase, istep,1);
            end
            
            if Add_data(1).ised == 1
                Add_data(array_index).Csed(:,count) = read_binary_file_pcui_profile(working_folder, fname_Sed, istep,params,0,0);
                Add_data(array_index).vCsed(:,count) = read_binary_file_pcui_profile(working_folder, fname_vCSed, istep,params,0,0);
                Add_data(array_index).CsedTotal1(1,count) = read_binary_file_pcui_value(working_folder,   fname_sedTotal1, istep,1);
                Add_data(array_index).CsedTotal2(1,count) = read_binary_file_pcui_value(working_folder,   fname_sedTotal2, istep,1);
                Add_data(array_index).BruntN(:,count) = read_binary_file_pcui_profile(working_folder, fname_BruntN, istep,params,0,0);
                Add_data(array_index).rhoMean(:,count) = read_binary_file_pcui_profile(working_folder, fname_rhoMean, istep,params,0,0);
                Add_data(array_index).rhoSqrMean(:,count) = read_binary_file_pcui_profile(working_folder, fname_rhoSqrMean, istep,params,0,0);
                Add_data(array_index).rhoPrimeMean(:,count) = read_binary_file_pcui_profile(working_folder, fname_rhoPrimeMean, istep,params,0,0);
            end
            
            if Add_data(1).ieddy == 1
                Add_data(array_index).vst(:,count) = read_binary_file_pcui_profile(working_folder, fname_vstMean, istep,params,0,0);
                Add_data(array_index).rr(:,count) = read_binary_file_pcui_profile(working_folder, fname_rruMean, istep,params,0,0);
            end
            count =count+1;
        end
    end
end

%%%%%%%Velocity Plotting%%%%%%%%%%%%%%%%%
fig1 = figure;
set(0,'DefaultAxesFontSize',16)
set(fig1, 'Position', [100, 100, 1049, 895]);
[mtemp,ntemp]=size(Add_data(1).uDepthAve);
plot((1:ntemp)*iskip-iskip-1,Add_data(array_index).uDepthAve(1,1:end))
xlab = xlabel('\textbf{$n$}');
set(xlab,'interpreter','Latex','FontSize',16,'FontWeight','bold')
ylab = ylabel('\textbf{$\overline{u}$}');
set(ylab,'interpreter','Latex','FontSize',16,'FontWeight','bold')

    
load('./DataFolder/dataStrat.mat','dataWC','params')
dataWC(replace_profile).dataFolder = allFolder;
dataWC(replace_profile).umean = Add_data(array_index).umean(:,:);
dataWC(replace_profile).vmean  = Add_data(array_index).vmean(:,:);
dataWC(replace_profile).wmean = Add_data(array_index).wmean(:,:);
dataWC(replace_profile).uTurb = Add_data(array_index).uTurb(:,:);
dataWC(replace_profile).vTurb = Add_data(array_index).vTurb(:,:);
dataWC(replace_profile).wTurb = Add_data(array_index).wTurb(:,:);
dataWC(replace_profile).uvRey = Add_data(array_index).uvRey(:,:);
dataWC(replace_profile).uwRey = Add_data(array_index).uwRey(:,:);
dataWC(replace_profile).vwRey = Add_data(array_index).vwRey(:,:);
dataWC(replace_profile).dissipation = Add_data(array_index).dissipation(:,:);
dataWC(replace_profile).uDepthAve = Add_data(array_index).uDepthAve(1,:);
dataWC(replace_profile).time = Add_data(array_index).time(1,:);
dataWC(replace_profile).kount = Add_data(array_index).kount(1,:);
dataWC(replace_profile).Sdrive = Add_data(array_index).Sdrive(1,:);
dataWC(replace_profile).totalKinetic = Add_data(array_index).totalKinetic(1,:);
%dataCurrentOnly(replace_profile).cfl = Add_data(array_index).cfl(1,:);
dataWC(replace_profile).dt = Add_data(1).dt;
dataWC(replace_profile).nsave = Add_data(1).nsave;
dataWC(replace_profile).nsavePro = Add_data(1).nsavePro;
dataWC(replace_profile).dpdxSteady = Add_data(1).dpdxSteady;
dataWC(replace_profile).dpdxWave = Add_data(1).dpdxWave;
dataWC(replace_profile).iTKE = Add_data(1).iTKE;
dataWC(replace_profile).waves = Add_data(1).waves;
dataWC(replace_profile).Twave = Add_data(1).Twave;
dataWC(replace_profile).ised = Add_data(1).ised;

if Add_data(1).ised == 1
    dataWC(replace_profile).rhoSed = Add_data(1).rhoSed;
    dataWC(replace_profile).tauCrit = Add_data(1).tauCrit;
    dataWC(replace_profile).Ased = Add_data(1).Ased;
    dataWC(replace_profile).nSed = Add_data(1).nsed;
    dataWC(replace_profile).ws = Add_data(1).ws;
    dataWC(replace_profile).ak = Add_data(1).ak;
    dataWC(replace_profile).Csed = Add_data(array_index).Csed(:,:);
    dataWC(replace_profile).vCsed = Add_data(array_index).vCsed(:,:);
    dataWC(replace_profile).CsedTotal1 = Add_data(array_index).CsedTotal1(1,:);
    dataWC(replace_profile).CsedTotal2 = Add_data(array_index).CsedTotal2(1,:);
    dataWC(replace_profile).BruntN = Add_data(array_index).BruntN(:,:);
    dataWC(replace_profile).rhoMean = Add_data(array_index).rhoMean(:,:);
    dataWC(replace_profile).rhoSqrMean = Add_data(array_index).rhoSqrMean(:,:);
    dataWC(replace_profile).rhoPrimeMean = Add_data(array_index).rhoPrimeMean(:,:);
end

if Add_data(1).ieddy == 1
    dataWC(replace_profile).vst = Add_data(array_index).vst(:,:);
    dataWC(replace_profile).rr = Add_data(array_index).rr(:,:);
end
if Add_data(1).iTKE == 1
    dataWC(replace_profile).proDepth = Add_data(array_index).proDepth(1,:);
    dataWC(replace_profile).disDepth = Add_data(array_index).disDepth(1,:);
end
if Add_data(1).waves == 1
    dataWC(replace_profile).phase = Add_data(array_index).phase(1,:);
end


save('./DataFolder/dataStrat.mat','dataWC','params','-v7.3')



