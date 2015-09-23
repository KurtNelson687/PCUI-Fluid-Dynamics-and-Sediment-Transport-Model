function [S2] = calculate_binary_S2(u,metrics_2d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates vertical shear using binary output from 
% PCUI. x,y are the grid points (cell centers) and u is the laterally (z)
% averaged density at the cell centers, both including 1 halo cell. 
%
% Bobby Arthur
% May 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate density gradient drho_dy = drho/dx_2
du_dy = 1./metrics_2d.J.* ...
        ( metrics_2d.XI_Y.*( u(3:end,2:end-1) - u(1:end-2,2:end-1) )/2 ...
        + metrics_2d.ET_Y.*( u(2:end-1,3:end) - u(2:end-1,1:end-2) )/2 );

%Calculate N^2
S2 = du_dy.^2;

end
