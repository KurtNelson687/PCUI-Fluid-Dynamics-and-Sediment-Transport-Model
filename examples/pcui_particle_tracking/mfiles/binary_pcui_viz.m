% Displays PCUI simulation data from binary output files
clear all; clc; close all;

% PLOTTING OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timestep = 0; %timestep to plot
delta_ts = 0; %averaging
FIGURE_ON = 1; %figure visible?
    print_ext = '-dpng'; %image file type
    print_res = '-r200'; %image resolution

display_grid = 0;
display_density = 0;
display_velocity = 0;
display_scalar = 0;
display_pressure = 0;
display_vst = 0;
display_akst = 0;
display_diss_sgs = 0;
display_density_isosurface = 0;
    rho_iso = 1;
display_omega_1_isosurface = 0;
    omega_1_iso = 1;

x_loc = 0; %if 0, west boundary
y_loc = 0; %if 0, centerline

plot_xz = 1; %x-z plot
plot_yz = 0; %y-z plot

plot_quiver = 0;
show_grid_lines = 0; %false = shading flat

working_folder = '../';
save_folder =    '../figs';
% These are the output files with the filenames stripped out of extensions 
% (extensions are chosen automatically based on number of processors).
fname_xyz = 'xyz';
fname_rho = 'output_S';
fname_phi = 'output_phi';
fname_uvw = 'output_UVW';
filename_xpart = 'output_xPart.dat';
filename_upart = 'output_uPart.dat';
fname_vst = 'output_vst_o';
fname_akst = 'output_akst_o';
fname_diss_sgs = 'output_diss_sgs';

% OTHER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the file containing the parameter definition
ftext = fileread(fullfile(working_folder, 'io.f'));
params.dt = variable_value_pcui('dtime',ftext);
params.molecular_viscosity = variable_value_pcui('vis',ftext);
params.eddy_viscosity = variable_value_pcui('ak',ftext);
params.nsteps = variable_value_pcui('nstep',ftext);
params.nsave = variable_value_pcui('nsave',ftext);

% read the file containing the domain definition
ftext = fileread(fullfile(working_folder, 'cavity.f'));
params.bx = variable_value_pcui('bx',ftext); x_length = params.bx;
params.by = variable_value_pcui('by',ftext); z_length = params.by;
params.bz = variable_value_pcui('bz',ftext); y_length = params.bz;

% read the file containing the grid size and processor definitions
ftext = fileread(fullfile(working_folder, 'size.inc'));
params.ni = variable_value_pcui('ni',ftext); Nx = params.ni;
params.nj = variable_value_pcui('nj',ftext); Nz = params.nj;
params.nk = variable_value_pcui('nk',ftext); Ny = params.nk;
params.px = variable_value_pcui('px',ftext);
params.py = variable_value_pcui('py',ftext);
params.pz = variable_value_pcui('pz',ftext);

%Find correct istep value
n = [0, params.nsave:params.nsave:params.nsteps, params.nsteps+1];
istep = find(n==timestep,1);

% GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the files containing the grid definition and assemble into array
% (includes 2 halo cells for each proc)
[xh,yh,zh] = read_binary_file_pcui(working_folder, fname_xyz, 1, params,1,1);

%extract physical grid points including 1 halo cell
nni = params.ni/params.px+4;
nnj = params.nj/params.py+4;
nnk = params.nk/params.pz+4;
exi = []; exj = []; exk = [];

if params.px == 1
    exi = 2:nni-1;
else
    for i=0:params.px-1
        if i==0
            exi = [exi,2:nni-2];
        elseif i==params.px-1
            exi = [exi,i*nni+3:(i+1)*nni-1];
        else
            exi = [exi,i*nni+3:(i+1)*nni-2];
        end
    end
end

if params.py == 1
    exj = 2:nnj-1;
else
    for j=0:params.py-1
        if j==0
            exj = [exj,2:nnj-2];
        elseif j==params.py-1
            exj = [exj,j*nnj+3:(j+1)*nnj-1];
        else
            exj = [exj,j*nnj+3:(j+1)*nnj-2];
        end
    end
end

if params.pz == 1
    exk = 2:nnk-1;
else
    for k=0:params.pz-1
        if k==0
            exk = [exk,2:nnk-2];
        elseif k==params.pz-1
            exk = [exk,k*nnk+3:(k+1)*nnk-1];
        else
            exk = [exk,k*nnk+3:(k+1)*nnk-2];
        end
    end
end

%extract values and rearrange
xh = xh(exi,exj,exk); xh = permute(xh,[1 3 2]);
temp = yh(exi,exj,exk);
yh = zh(exi,exj,exk); yh = permute(yh,[1 3 2]);
zh = temp; zh = permute(zh,[1 3 2]);

%without halo
x = xh(2:end-1,2:end-1,2:end-1);
y = yh(2:end-1,2:end-1,2:end-1);
z = zh(2:end-1,2:end-1,2:end-1);

if(x_loc)
    x_slice = find(squeeze(x(:,1,1))>=x_loc,1,'first');
else
    x_slice = 1;
end

if(y_loc)
    y_slice = find(squeeze(y(1,:,1))>=y_loc,1,'first');
else
    y_slice = round(Ny/2); 
end

if(plot_xz)
   x_xz = squeeze(x(:,y_slice,:)); 
   z_xz = squeeze(z(:,y_slice,:)); 
end

if(plot_yz)
    y_yz = squeeze(y(x_slice,:,:));
    z_yz = squeeze(z(x_slice,:,:));
end

if(display_grid)  
    if(plot_xz)
        grid_fig_xz = figure;
        hold all;
        if(~FIGURE_ON)
            set(grid_fig_xz,'visible','off');
        end
        plot(x_xz(:,:), z_xz(:,:), 'k.'); 
        axis equal;
        axis([0 x_length -z_length 0]);
        xlabel('x [m]');
        ylabel('z [m]');
        title('Computational grid, x-z slice');
        %print(grid_fig_xz,print_ext,print_res,'grid_xz');
        %saveas(grid_fig_xz,'grid_xz.fig');
    end
    
    if(plot_yz)
        grid_fig_yz = figure;
        hold all;
        if(~FIGURE_ON)
            set(grid_fig_yz,'visible','off');
        end
        plot(y_yz(:,:), z_yz(:,:), 'k.');         
        axis equal;
        axis([0 y_length -z_length 0]);
        xlabel('y [m]');
        ylabel('z [m]');
        title('Computational grid, y-z slice');
        %print(grid_fig_yz,print_ext,print_res,'grid_yz');
        %saveas(grid_fig_yz,'grid_yz.fig');
    end
end

% DENSITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_density)
    rho = read_binary_file_pcui(working_folder, fname_rho, istep, ...
                                 params, 0,0);     
    rho = permute(rho,[1 3 2]);
   
   if(plot_xz)
       rho_xz = squeeze(rho(:,y_slice,:));
       rho_fig_xz = figure;
       hold all;
       if(~FIGURE_ON)
           set(rho_fig_xz,'visible','off');
       end
       pcolor(x_xz,z_xz,rho_xz);
       colorbar;
%        contour(x_xz,z_xz,rho_xz,[1, 1],'k');
       axis equal;
%        axis([0 x_length -z_length 0]);
       xlabel('x [m]');
       ylabel('z [m]');
       title('Density \Delta\rho/\rho_0, x-z slice');
       if(~show_grid_lines)
           shading flat;
       end
       if (plot_yz)
           hold on;
           plot([x_loc x_loc],[z(x_slice,1,1) z(x_slice,1,end)],'k-','LineWidth',2);
           hold off;
       end
       %print(rho_fig_xz,print_ext,print_res,'density_xz');
       %saveas(rho_fig_xz,'density_xz.png');
   end

   if(plot_yz)
       rho_yz = squeeze(rho(x_slice,:,:));
       rho_fig_yz = figure;
       hold all;
       if(~FIGURE_ON)
           set(rho_fig_yz,'visible','off');
       end
       pcolor(y_yz,z_yz,rho_yz);
       colorbar;
       axis image;
%        axis([0 y_length -z_length 0]);
       xlabel('y [m]');
       ylabel('z [m]');
       title('Density \Delta\rho/\rho_0, y-z slice');
       if(~show_grid_lines)
           shading flat;
       end
       %print(rho_fig_yz,print_ext,print_res,'density_yz');
       %saveas(rho_fig_yz,'density_yz.fig');
   end
end

% VELOCITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_velocity)
    [u,v,w] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
                                      params, 1,0);     
    u = permute(u,[1 3 2]);
    temp = v;
    v = w; v = permute(v,[1 3 2]);
    w = temp; w = permute(w,[1 3 2]);
    
    if(plot_xz)
        u_xz = squeeze(u(:,y_slice,:));
        v_xz = squeeze(v(:,y_slice,:));
        w_xz = squeeze(w(:,y_slice,:));
        velocity_fig_xz = figure;
        hold all;
        if(~FIGURE_ON)
           set(velocity_fig_xz,'visible','off');
        end
        pcolor(x_xz,z_xz,u_xz);
        colorbar;
        if(plot_quiver)
            quiver(x_xz,z_xz,u_xz,w_xz,'k');
        end
        axis equal;
%         axis([0 x_length -z_length 0]);
        xlabel('x [m]');
        ylabel('z [m]');
        title('Velocity u, x-z slice');
        if(~show_grid_lines)
            shading flat;
        end
        %print(velocity_fig_xz,print_ext,print_res,'velocity_xz');
        %saveas(velocity_fig_xz,'velocity_xz.fig');
    end
    
    if(plot_yz)
        u_yz = squeeze(u(x_slice,:,:));
        v_yz = squeeze(v(x_slice,:,:));
        w_yz = squeeze(w(x_slice,:,:));
        velocity_fig_yz = figure;
        hold all;
        if(~FIGURE_ON)
            set(velocity_fig_yz,'visible','off');
        end
        pcolor(y_yz,z_yz,u_yz);
        colorbar;
        if(plot_quiver)
            quiver(y_yz,z_yz,v_yz,w_yz);
        end
        axis equal;
%         axis([0 y_length -z_length 0]);
        xlabel('y [m]');
        ylabel('z [m]');
        title('Velocity u, y-z slice');
        if(~show_grid_lines)
            shading flat;
        end
        %print(velocity_fig_yz,print_ext,print_res,'velocity_yz');
        %saveas(velocity_fig_yz,'velocity_yz.fig');
    end  
end

% SCALAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_scalar)
    phi = read_binary_file_pcui(working_folder, fname_phi, istep, ...
                                 params, 0,0);     
    phi = permute(phi,[1 3 2]);

   if(plot_xz)
       phi_xz = squeeze(phi(:,y_slice,:));
       phi_fig_xz = figure;
       hold all;
       if(~FIGURE_ON)
          set(phi_fig_xz,'visible','off');
       end
       pcolor(x_xz,z_xz,phi_xz);
       colorbar;
%        caxis([0 1]);
       if(plot_quiver)
           quiver(x_xz,z_xz,u_xz,w_xz);
       end
       axis equal;
       axis([0 x_length -z_length 0]);
       xlabel('x [m]');
       ylabel('z [m]');
       title('Scalar \phi, x-z slice');
       if(~show_grid_lines)
           shading flat;
       end
       %print(phi_fig_xz,print_ext,print_res,'scalar_xz');
       %saveas(phi_fig_xz,'scalar_xz.fig');
   end
   
   if(plot_yz)
       phi_yz = squeeze(phi(x_slice,:,:));
       phi_fig_yz = figure;
       hold all;
       if(~FIGURE_ON)
          set(phi_fig_yz,'visible','off');
       end
       pcolor(y_yz,z_yz,phi_yz);
       colorbar;
%        caxis([0 1]);
       if(plot_quiver)
           quiver(y_yz,z_yz,v_yz,w_yz);
       end
       axis equal;
       axis([0 y_length -z_length 0]);
       xlabel('y [m]');
       ylabel('z [m]');
       title('Scalar \phi, y-z slice');
       if(~show_grid_lines)
           shading flat;
       end
       %print(phi_fig_yz,print_ext,print_res,'scalar_yz');
       %saveas(phi_fig_yz,'scalar_yz.fig');
   end
end

% EDDY VISCOSITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_vst)
    vst = read_binary_file_pcui(working_folder, fname_vst, istep, ...
                                 params, 0,0);     
    vst = permute(vst,[1 3 2]);
   
   if(plot_xz)
       vst_xz = squeeze(vst(:,y_slice,:));
       vst_fig_xz = figure;
       hold all;
       if(~FIGURE_ON)
           set(vst_fig_xz,'visible','off');
       end
       pcolor(x_xz,z_xz,vst_xz);
       colorbar;
%        contour(x_xz,z_xz,rho_xz,[1, 1],'k');
       axis equal;
%        axis([0 x_length -z_length 0]);
       xlabel('x [m]');
       ylabel('z [m]');
       title('Eddy Viscosity, x-z slice');
       if(~show_grid_lines)
           shading flat;
       end
       if (plot_yz)
           hold on;
           plot([x_loc x_loc],[z(x_slice,1,1) z(x_slice,1,end)],'k-','LineWidth',2);
           hold off;
       end
       %print(rho_fig_xz,print_ext,print_res,'density_xz');
       %saveas(rho_fig_xz,'density_xz.png');
   end
end

% EDDY DIFFUSIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_akst)
    akst = read_binary_file_pcui(working_folder, fname_akst, istep, ...
                                 params, 0,0);     
    akst = permute(akst,[1 3 2]);
   
   if(plot_xz)
       akst_xz = squeeze(akst(:,y_slice,:));
       akst_fig_xz = figure;
       hold all;
       if(~FIGURE_ON)
           set(akst_fig_xz,'visible','off');
       end
       pcolor(x_xz,z_xz,akst_xz);
       colorbar;
%        contour(x_xz,z_xz,rho_xz,[1, 1],'k');
       axis equal;
%        axis([0 x_length -z_length 0]);
       xlabel('x [m]');
       ylabel('z [m]');
       title('Eddy Diffusivity, x-z slice');
       if(~show_grid_lines)
           shading flat;
       end
       if (plot_yz)
           hold on;
           plot([x_loc x_loc],[z(x_slice,1,1) z(x_slice,1,end)],'k-','LineWidth',2);
           hold off;
       end
       %print(rho_fig_xz,print_ext,print_res,'density_xz');
       %saveas(rho_fig_xz,'density_xz.png');
   end
end

% SGS DISSIPATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_diss_sgs)
    diss_sgs = read_binary_file_pcui(working_folder, fname_diss_sgs, istep, ...
                                 params, 0,0);     
    diss_sgs = permute(diss_sgs,[1 3 2]);
   
   if(plot_xz)
       diss_sgs_xz = squeeze(diss_sgs(:,y_slice,:));
       diss_sgs_fig_xz = figure;
       hold all;
       if(~FIGURE_ON)
           set(diss_sgs_fig_xz,'visible','off');
       end
       pcolor(x_xz,z_xz,diss_sgs_xz);
       colorbar;
%        contour(x_xz,z_xz,rho_xz,[1, 1],'k');
       axis equal;
%        axis([0 x_length -z_length 0]);
       xlabel('x [m]');
       ylabel('z [m]');
       title('SGS Dissipation, x-z slice');
       if(~show_grid_lines)
           shading flat;
       end
       if (plot_yz)
           hold on;
           plot([x_loc x_loc],[z(x_slice,1,1) z(x_slice,1,end)],'k-','LineWidth',2);
           hold off;
       end
       %print(rho_fig_xz,print_ext,print_res,'density_xz');
       %saveas(rho_fig_xz,'density_xz.png');
   end
end

% PRESSURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_pressure)    
    [p] = load_binary_pressure(directory, timestep,delta_ts, npx,npy,npz, l_ni_h1,l_nj_h1,l_nk_h1);

    if(plot_xz)
       p_xz = squeeze(p(:,y_slice,:));
       p_fig_xz = figure;
       hold all;
       if(~FIGURE_ON);
           set(p_fig_xz,'visible','off');
       end
       pcolor(x_xz,z_xz,p_xz);
       colorbar;
       axis equal;
       axis([0 x_length -z_length 0]);
       xlabel('x [m]');
       ylabel('z [m]');
       title('Pressure, x-z slice');
       if(~show_grid_lines)
           shading flat;
       end
       %print(p_fig_xz,print_ext,print_res,'pressure_xz');
       %saveas(p_fig_xz,'pressure_xz.fig');
   end
   
   if(plot_yz)
       phi_yz = squeeze(p(x_slice,:,:));
       p_fig_yz = figure;
       hold all;
       if(~FIGURE_ON);
           set(p_fig_yz,'visible','off');
       end
       pcolor(y_yz,z_yz,p_yz);
       colorbar;
       axis equal;
       axis([0 y_length -z_length 0]);
       xlabel('y [m]');
       ylabel('z [m]');
       title('Pressure, y-z slice');
       if(~show_grid_lines)
           shading flat;
       end
       %print(p_fig_yz,print_ext,print_res,'pressure_yz');
       %saveas(p_fig_yz,'pressure_yz.fig');
   end
end

% DENSITY ISOSURFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(display_density_isosurface)
    z_slice = round(Nz/2);
    xz = squeeze(x(:,:,z_slice)); yz = squeeze(y(:,:,z_slice)); 
    depth = abs(squeeze(z(:,:,1)));

    rho = read_binary_file_pcui(working_folder, fname_rho, istep, ...
                                 params, 0,0);     
    rho = permute(rho,[1 3 2]);

    isoplot = figure;
    hold all;
    if(~FIGURE_ON)
        set(isoplot,'visible','off');
    end
    surf(xz,yz,-depth-.005,'FaceColor','k','EdgeColor','none'); 
    pat = patch(isosurface(x,y,z,rho,rho_iso,rho));
    set(pat,'FaceColor','interp','EdgeColor','none');
    set(pat,'FaceColor','r');
    
    if(display_omega_1_isosurface)
        [u,v,w] = read_binary_file_pcui(working_folder, fname_uvw, istep, ...
                                          params, 1,1);     
        u = u(exi,exj,exk); u = permute(u,[1 3 2]);
        temp = v(exi,exj,exk);
        v = w(exi,exj,exk); v = permute(v,[1 3 2]);
        w = temp; w = permute(w,[1 3 2]);
        metrics = calculate_binary_metrics(xh,yh,zh);
        omega_1 = calculate_binary_vorticity(metrics,v,w);
        pat_pos = patch(isosurface( ...
                     x, ...
                     y, ...
                     z, ...
                     omega_1,omega_1_iso,omega_1) );
        set(pat_pos,'FaceColor','b','EdgeColor','none');
        pat_neg = patch(isosurface( ...
                     x, ...
                     y, ...
                     z, ...
                     omega_1,-omega_1_iso,omega_1) );
        set(pat_neg,'FaceColor','g','EdgeColor','none');
    end
    
    daspect([1,1,1])
    view([0.04, -0.04, 0.04]); axis tight
    camlight 
    lighting gouraud
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    title(['Isosurface of \Delta\rho/\rho_0=',num2str(rho_iso)]);
    %print(isoplot,print_ext,print_res,'isoplot_3D');
    %saveas(isoplot,'isoplot_3D.fig');
end

% PARTICLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xpart = read_binary_particles_pcui(working_folder, filename_xpart, istep, params);
plot3(xpart(:,1),xpart(:,3),xpart(:,2),'k.');
% plot(xpart(:,1),xpart(:,2),'k.');
% plot(xpart(:,1),xpart(:,3),'k.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('Complete');
