function gs = calculate_binary_gradient_squared_2d(rho,metrics_2d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the squared gradient of rho in two dimensions (x,y)
%
% | grad(rho) |^2 = d(rho)/dx_j*(drho)/dx_j
%
% Bobby Arthur
% May 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate components of density gradient s = d(rhop)/dx_j
d_1 = 1./metrics_2d.J.* ...
     ( metrics_2d.XI_X.*( rho(3:end,2:end-1) - rho(1:end-2,2:end-1) )/2 ...
     + metrics_2d.ET_X.*( rho(2:end-1,3:end) - rho(2:end-1,1:end-2) )/2 );
 
d_2 = 1./metrics_2d.J.* ...
     ( metrics_2d.XI_Y.*( rho(3:end,2:end-1) - rho(1:end-2,2:end-1) )/2 ...
     + metrics_2d.ET_Y.*( rho(2:end-1,3:end) - rho(2:end-1,1:end-2) )/2 );

gs = d_1.^2 + d_2.^2;

end