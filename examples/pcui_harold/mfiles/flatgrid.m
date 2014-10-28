clear all; close all; clc;
load '/home/barthur/zang/grids/grid_1152x128x16_r102_w1_zstretch_s218_f128.mat';
xnew = repmat(x(1:612,1,130),[1 20 260]);
ynew = repmat(y(3,:,130),[612 1 260]);

%znew
Nz = 256; %128
Hint = 0.26; %Hint-H=-0.3
H = 0.56;
dz = H/Nz;
% znew = [0:dz:H];
znew = -H-1.5*dz:dz:1.5*dz;
znew = reshape(znew,[1 1 length(znew)]);
znew = repmat(znew,[612 20 1]);

x=xnew; y=ynew; z=znew;
% save('2ch_flat_part.mat','x','y','z');

% plot3(xnew(:),ynew(:),znew(:));

% % Height of resolved region
% Hr = 0.35;
% 
% % Number of grid points in resolved region
% Nzr = 100;
% dzr = Hr/Nzr;
% 
% % Maximum acceptable stretching factor. The code will exit if the
% % stretching factor needed to give the refined region exceeds rmax.
% rmax = 1.1;
% 
% %Stretch grid BOTTOM
% Hs = H-Hr;
% alpha = Hs*(Nzr/Hr);
% Ns = Nz - Nzr;
% 
% fsolveopts = optimset('Diagnostics','off');
% r = fsolve(@(r,Ns,alpha) r^Ns - alpha*r - 1 + alpha,1.1,...
% 	   fsolveopts,Ns,alpha);
% if(r<1 | r>rmax)
%   error(sprintf('Stretching factor is r=%.3f not in range [1,%.2f]',r,rmax));
% else
%   fprintf('Stretching factor is %.2f\n',r);
% end
% 
% zps = zeros(1,Nz+1);
% zps(1:Nzr+1) = [0:dzr:Hr]';
% 
% dz=dzr;
% for j=Nzr+1:Nz
%   zps(j+1) = (zps(j)+dz);
%   dz=r*dz;
% end
% 
% zps = zps - H;
% 
% %plot
% figure;
% plot(zeros(size(zps)),zps,'k.');
% figure;
% plot(diff(zps),zps(1:end-1),'k-');
% 
% %Finish grid
% zps = [zps(1)-2*diff(zps(1:2)), zps(1)-diff(zps(1:2)), ...
%        zps, ...
%        zps(end)+diff(zps(end-1:end)), zps(end)+2*diff(zps(end-1:end))];
% zps = resize(zps,[1 1 length(zps)]);
% znew = repmat(zps,[1156 20 1]);

% %Stretch grid MIDDLE
% Hs = (H-Hr)/2;
% alpha = Hs*(Nzr/Hr);
% Ns = (Nz - Nzr)/2;
% 
% fsolveopts = optimset('Diagnostics','off');
% r = fsolve(@(r,Ns,alpha) r^Ns - alpha*r - 1 + alpha,1.1,...
% 	   fsolveopts,Ns,alpha);
% if(r<1 | r>rmax)
%   error(sprintf('Stretching factor is r=%.3f not in range [1,%.2f]',r,rmax));
% else
%   fprintf('Stretching factor is %.2f\n',r);
% end
% 
% zps = zeros(1,Nz+1);
% zps(Ns+1:Ns+Nzr+1) = H/2+[-Hr/2:dzr:Hr/2]';
% 
% dz=dzr;
% for j=Ns+1:-1:2
%   zps(j-1) = (zps(j)-dz);
%   dz=r*dz;
% end
% 
% dz=dzr;
% for j=Ns+Nzr+1:Nz
%   zps(j+1) = (zps(j)+dz);
%   dz=r*dz;
% end
% 
% zps = zps - H;
% 
% %plot
% figure;
% plot(zeros(size(zps)),zps,'k.');
% figure;
% plot(diff(zps),zps(1:end-1),'k-');

% %Stretch grid INTERFACE
% %Bottom
% Hs = Hint-Hr/2;
% alpha = Hs*(Nzr/Hr);
% Ns = round((Nz - Nzr/2)*(Hint-Hr/2)/H);
% Ns_old = Ns;
% Hs_old = Hs;
% 
% fsolveopts = optimset('Diagnostics','off');
% r = fsolve(@(r,Ns,alpha) r^Ns - alpha*r - 1 + alpha,1.1,...
% 	   fsolveopts,Ns,alpha);
% if(r<1 | r>rmax)
%   error(sprintf('Stretching factor is r=%.3f not in range [1,%.2f]',r,rmax));
% else
%   fprintf('Stretching factor is %.2f\n',r);
% end
% 
% zps = zeros(1,Nz+1);
% zps(Ns+1:Ns+Nzr+1) = Hint+[-Hr/2:dzr:Hr/2]';
% 
% dz=dzr;
% for j=Ns+1:-1:2
%   zps(j-1) = (zps(j)-dz);
%   dz=r*dz;
% end
% 
% %Top
% Hs = H-Hint-Hr/2;
% alpha = Hs*(Nzr/Hr);
% Ns = Nz - Nzr - Ns_old;
% 
% fsolveopts = optimset('Diagnostics','off');
% r = fsolve(@(r,Ns,alpha) r^Ns - alpha*r - 1 + alpha,1.1,...
% 	   fsolveopts,Ns,alpha);
% if(r<1 | r>rmax)
%   error(sprintf('Stretching factor is r=%.3f not in range [1,%.2f]',r,rmax));
% else
%   fprintf('Stretching factor is %.2f\n',r);
% end
% 
% dz=dzr;
% for j=Ns+Nzr+1:Nz
%   zps(j+1) = (zps(j)+dz);
%   dz=r*dz;
% end
% 
% zps = zps - H;
% 
% %plot
% figure;
% plot(zeros(size(zps)),zps,'k.');
% figure;
% plot(diff(zps),zps(1:end-1),'k-');

