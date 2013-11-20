% Loads CNS simulation data from binary output files to .mat file
clear all; clc; close all;

working_folder = '/usr/var/tmp/barthur/';

% LOAD OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'midslice.mat';
timestep_initial = 0;
timestep_final = 10000; 

load_grid = 1;
load_density = 1;
load_velocity = 0;
load_scalar = 0;
load_pressure = 0;
load_potential_energy = 0;

lateral_average = 0;
mid_slice = 1;
    y_slice = 16;

% LOAD STATIC PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Choose variable to save
save_vars = '''Nx'',''Ny'',''Nz'',''params'',''t'',''n''';
if(load_grid)
    save_vars = [save_vars,',''x'',''y'',''z'''];
end
if(load_density)
    save_vars = [save_vars,',''rho'''];
end
if(load_velocity)
    save_vars = [save_vars,',''u'',''v'',''w'''];
end
if(load_scalar)
    save_vars = [save_vars,',''phi'''];
end
if(load_pressure)
    save_vars = [save_vars,',''p'''];
end
if(load_potential_energy)
    save_vars = [save_vars,',''Eb'',''Ep'',''Ek'',''epsilon_tot'''];
end
if(lateral_average)
    save_vars = [save_vars,',''lateral_average'''];
end

% read the file containing the parameter definition
ftext = fileread(fullfile(working_folder, 'io.f'));
params.dt = variable_value_pcui('dtime',ftext);
params.molecular_viscosity = variable_value_pcui('vis',ftext);
params.eddy_viscosity = variable_value_pcui('ak',ftext);
params.nsteps = variable_value_pcui('nstep',ftext);
params.nsave = variable_value_pcui('nsave',ftext);

% read the file containing the domain definition
ftext = fileread(fullfile(working_folder, 'cavity.f'));
params.bx = variable_value_pcui('bx',ftext); 
params.by = variable_value_pcui('by',ftext); 
params.bz = variable_value_pcui('bz',ftext); 

% read the file containing the grid size and processor definitions
ftext = fileread(fullfile(working_folder, 'size.inc'));
params.ni = variable_value_pcui('ni',ftext); Nx = params.ni;
params.nj = variable_value_pcui('nj',ftext); Nz = params.nj;
params.nk = variable_value_pcui('nk',ftext); Ny = params.nk;
params.px = variable_value_pcui('px',ftext);
params.py = variable_value_pcui('py',ftext);
params.pz = variable_value_pcui('pz',ftext);

% These are the output files with the filenames stripped out of extensions 
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
fname_rho = 'output_S';
fname_uvw = 'output_UVW';

%Grid
if(load_grid)
    [x,y,z] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,0);
    x = permute(x,[1 3 2]); 
    temp = permute(y,[1 3 2]);
    y = permute(z,[1 3 2]);
    z = temp;
    
    if(lateral_average || mid_slice)
        x = x(:,1,:);
        y = y(:,1,:);
        z = z(:,1,:);
    end
end

% LOAD TIME VARIABLE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Time
n = [0, params.nsave:params.nsave:params.nsteps, params.nsteps+1];
t = params.dt*n;

istart = find(n==timestep_initial,1);
iend = find(n==timestep_final,1);

%Initialize storage variables
if(lateral_average || mid_slice)
    rho = zeros(Nx,Nz,iend);
    u = rho; v = rho; w = rho; phi = rho; p = rho;
else
    rho = zeros(Nx,Ny,Nz,iend);
    u = zeros(Nx,Ny,Nz,iend); v = u; w = u; phi = rho; p = rho;
end

for istep=istart:iend
    str = sprintf(['istep = ',num2str(istep),', t = ',num2str(t(istep))]); 
    disp(str);
    
    %Density
    if(load_density)
        rho_n = read_binary_file_pcui(working_folder, fname_rho, istep, ...
                                     params, 0,0);     
        rho_n = permute(rho_n,[1 3 2]);
        
        if(lateral_average)
            rho(:,:,istep) = mean(rho_n,2);
        elseif(mid_slice)
            rho(:,:,istep) = rho_n(:,y_slice,:);
        else
            rho(:,:,:,istep) = rho_n;
        end
        clear 'rho_n';
    end
    
    %Velocity
    if(load_velocity)
        [u_n,v_n,w_n] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
                                      params, 1,0); 
        u_n = permute(u_n,[1 3 2]);
        temp = permute(v_n,[1 3 2]);
        v_n = permute(w_n,[1 3 2]);
        w_n = temp;
        
        if(lateral_average)
            u(:,:,istep) = mean(u_n,2);
            v(:,:,istep) = mean(v_n,2);
            w(:,:,istep) = mean(w_n,2);
        elseif(mid_slice)
            u(:,:,istep) = u_n(:,y_slice,:);
            v(:,:,istep) = v_n(:,y_slice,:);
            w(:,:,istep) = w_n(:,y_slice,:);
        else
            u(:,:,:,istep) = u_n;
            v(:,:,:,istep) = v_n;
            w(:,:,:,istep) = w_n;
        end
        clear 'u_n' 'v_n' 'w_n';
    end
    
%     %Scalar
%     if(load_scalar)
%         [phi_n] = load_binary_scalar(directory,tn,...
%             delta_ts/save_timestep_period,npx,npy,npz,l_ni_h2,l_nj_h2,l_nk_h2);
%         if(lateral_average)
%             phi(:,:,n) = mean(phi_n,2);
%         elseif(mid_slice)
%             phi(:,:,n) = phi_n(:,y_slice,:);
%         else
%             phi(:,:,:,n) = phi_n;
%         end
%         clear 'phi_n';
%     end
%     
%     %Pressure
%     if(load_pressure)    
%         [p_n] = load_binary_pressure(directory,tn,...
%             delta_ts, npx,npy,npz, l_ni_h1,l_nj_h1,l_nk_h1);
%         if(lateral_average)
%             p(:,:,n) = mean(p_n,2);
%         elseif(mid_slice)
%             p(:,:,n) = p_n(:,y_slice,:);
%         else
%             p(:,:,:,n) = p_n;
%         end
%         clear 'p_n';
%     end
end

%Squeeze out extra dimensions
t = t(istart:iend);
n = n(istart:iend);
rho = squeeze(rho);
u = squeeze(u);
v = squeeze(v);
w = squeeze(w);
phi = squeeze(phi);
p = squeeze(p);

%Save 
eval(['save(''',filename,''',',save_vars,')']);
display('Complete');

