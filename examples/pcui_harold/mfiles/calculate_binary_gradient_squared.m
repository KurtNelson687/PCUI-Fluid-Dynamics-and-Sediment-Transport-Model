function gs = calculate_binary_gradient_squared(rho,metrics)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the squared gradient of rho
%
% | grad(rho) |^2 = d(rho)/dx_j*(drho)/dx_j
%
% Bobby Arthur
% May 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate components of density gradient s = d(rhop)/dx_j
d_1 = 1./metrics.J.* ...
     ( metrics.XI_X.*( rho(3:end,2:end-1,2:end-1) - rho(1:end-2,2:end-1,2:end-1) )/2 ...
     + metrics.ET_X.*( rho(2:end-1,3:end,2:end-1) - rho(2:end-1,1:end-2,2:end-1) )/2 ...
     + metrics.ZT_X.*( rho(2:end-1,2:end-1,3:end) - rho(2:end-1,2:end-1,1:end-2) )/2 );
 
d_2 = 1./metrics.J.* ...
     ( metrics.XI_Y.*( rho(3:end,2:end-1,2:end-1) - rho(1:end-2,2:end-1,2:end-1) )/2 ...
     + metrics.ET_Y.*( rho(2:end-1,3:end,2:end-1) - rho(2:end-1,1:end-2,2:end-1) )/2 ...
     + metrics.ZT_Y.*( rho(2:end-1,2:end-1,3:end) - rho(2:end-1,2:end-1,1:end-2) )/2 );
 
d_3 = 1./metrics.J.* ...
     ( metrics.XI_Z.*( rho(3:end,2:end-1,2:end-1) - rho(1:end-2,2:end-1,2:end-1) )/2 ...
     + metrics.ET_Z.*( rho(2:end-1,3:end,2:end-1) - rho(2:end-1,1:end-2,2:end-1) )/2 ...
     + metrics.ZT_Z.*( rho(2:end-1,2:end-1,3:end) - rho(2:end-1,2:end-1,1:end-2) )/2 );

gs = d_1.^2 + d_2.^2 + d_3.^2;

end