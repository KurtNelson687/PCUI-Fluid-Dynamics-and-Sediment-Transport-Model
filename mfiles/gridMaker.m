%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kurt Nelson, 3/8/2016
% 
% This script makes a grid that is stretched in the vertical direction(y).
% I use this to examine my PCUI grids before running a simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all; clear all;clc

NJ=10*16;
NK=40*16;
NI=40*16;
max_aspect_ratio=25;
num_wall_unit_bed = 1;
stretch_factor = 5/100;
rho = 1000;
numPoints = NI*NJ*NK

nu = 1*10^-6; %kinematic viscosity
u_star = 0.01; %friction velocity
wall_unit = nu/u_star; %length of wall unit
dy_min = wall_unit*num_wall_unit_bed; %dy at bed
dx= 40*wall_unit;%dy_min*max_aspect_ratio; %dx (constant)
dz= 20*wall_unit;%dx;%dz (constant)

stop_stretch = floor(log(dx/dy_min)/log(1+stretch_factor)+1); %cell number where vertical stretching stops

%dy of each cell
for i= 1:NJ
    if i==1
        dy(i)=dy_min;
    elseif i>stop_stretch
        dy(i)=dx;
    else
        dy(i)=dy_min*(1+stretch_factor)^(i-1);
    end
    aspect_ratio(i)=dx/dy(i);
end

%y coordinate of cell center (note: y is vertical coordinate in PCUI)
y(1)=dy_min/2;
for i= 1:NJ-1
    growth_rate(i) = dy(i+1)/dy(i);
    y(i+1)=y(i)+(dy(i)+dy(i+1))/2;
end

%x coordinate of cell center (streamwise direction)
x(1)=dx/2;
for i=1:NI-1
    x(i+1)=x(i)+dx;
end

%z coordinate of cell center (cross-stream direction)
z(1)=dz/2;
for i=1:NK-1
    z(i+1)=z(i)+dz;
end

[X,Z,Y]=meshgrid(x,z,y);

% figure
% set(0,'DefaultAxesFontSize',16)
% xplot(:,:) = X(1,:,:);
% yplot(:,:) = Y(1,:,:);
% plot(xplot,yplot,'k',xplot',yplot','k')
% xlabel('x (m)')
% ylabel('y (m)')
% axis equal
% set(gca,'Xlim',[-max(x)/20 max(x)+max(x)/20])

domain_height=sum(dy)
domain_length_ratio=dx*NI/domain_height
domain_width_ratio=dz*NK/domain_height

dpdx = u_star^2*rho/(sum(dy)+0.5*dy(end))
%%% Modify method below if you want logarithmic stretching. I'm giving up
%%% on it for now. Its not as flexable as I thought.
% % r= [10^-5,10^-4,10^-3,10^-2,10^-1,10^-0]
% % value=depth-dx./log(r.*NJ).*log(factorial(NJ).*r.^NJ);
% % plot(r,value); 
% f=@(r) depth-dx/log(r*NJ)*log(factorial(NJ)*r^NJ); 
% 
% b=fzero(f,10^-3);
% figure
% plot(1:100,log(1:100))
% for i=1:NJ
%    dy(i)=dx*log(b*i)/log(b*NJ); 
% end
% 
% sum(dy)
% figure
% plot(1:NJ,dy)