clear all; close all; clc;
load gridtest.mat
xgrid = x;
zgrid = z;
load uvwrho_6000.mat
metrics = calculate_binary_metrics(x,y,z);
omega1 = calculate_binary_vorticity(metrics,v,w);

%%
nx=find(squeeze(x(:,1,end))>3.81,1); %3.951
nz=find(squeeze(z(nx,1,:))>-.2892,1); %-0.3901

figure;
pcolor(squeeze(y(nx,:,:)),squeeze(z(nx,:,:)),squeeze(w(nx,:,:))); 
caxis manual; colorbar; shading flat;
hold on;
plot(squeeze(y(nx,:,nz)),squeeze(z(nx,:,nz)),'k')
% contour(squeeze(y(nx,2:end-1,2:end-1)),squeeze(z(nx,2:end-1,2:end-1)),squeeze(omega1(nx,:,:)),10,'k');
% contour(squeeze(y(nx,:,:)),squeeze(z(nx,:,:)),squeeze(rho(nx,:,:)),10,'k');
axis image;
xlabel('y [m]');
ylabel('z [m]');
title('w(y,z) [m/s]')

figure;
plot(squeeze(y(nx,:,nz)),squeeze(w(nx,:,nz)),'k.-');

figure;
subplot(2,1,1)
plot(squeeze(y(nx,2:end-1,nz)),squeeze(omega1(nx,:,nz)),'k.-');
ylabel('\omega_1 [s^{-1}]')
subplot(2,1,2)
plot(squeeze(y(nx,:,nz)),squeeze(w(nx,:,nz)),'k.-');
xlabel('y [m]')
ylabel('w [m/s]')

figure;
hold on;
pcolor(squeeze(x(:,32,:)),squeeze(z(:,32,:)),squeeze(rho(:,32,:))); 
shading flat;
plot(squeeze(x(nx,1,:)),squeeze(z(nx,1,:)),'k-')
% plot(squeeze(xgrid(3:end-2,1,3:end-2)),squeeze(zgrid(3:end-2,1,3:end-2)),'k.')
axis image;
xlabel('x [m]')
ylabel('z [m]')

figure;
pcolor(squeeze(y(nx,2:end-1,2:end-1)),squeeze(z(nx,2:end-1,2:end-1)),squeeze(omega1(nx,:,:)));
axis equal;
colorbar;