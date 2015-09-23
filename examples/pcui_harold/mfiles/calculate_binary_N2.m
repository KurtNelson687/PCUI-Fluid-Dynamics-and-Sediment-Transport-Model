function [N2] = calculate_binary_N2(rho,g,metrics_2d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates buoyancy frequency using binary output from 
% PCUI. x,y are the grid points (cell centers) and rho is the laterally (z)
% averaged density at the cell centers, both including 1 halo cell. The 
% buoyancy frequency is calculated as 
%
% N^2 = g*(-drho/dy) (y is vertical)
%
% Bobby Arthur
% December 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate density gradient drho_dy = drho/dx_2
drho_dy = 1./metrics_2d.J.* ...
        ( metrics_2d.XI_Y.*( rho(3:end,2:end-1) - rho(1:end-2,2:end-1) )/2 ...
        + metrics_2d.ET_Y.*( rho(2:end-1,3:end) - rho(2:end-1,1:end-2) )/2 );

%Calculate N^2
N2 = -g*drho_dy;

end
