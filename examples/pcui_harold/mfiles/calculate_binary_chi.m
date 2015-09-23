function chi = calculate_binary_chi(rhop,kappa,metrics)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates laterally-averaged dissipation of scalar variance chi.
%
% chi = 2*kappa*d(rhop)/dx_j*(drhop)/dx_j
%
% where rhop is the density fluctuation about the laterval average.
%
% Bobby Arthur
% May 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate components of density gradient s = d(rhop)/dx_j
d_1 = 1/metrics.J.* ...
     ( metrics.XI_X.*( rhop(3:end,2:end-1,2:end-1) - rhop(1:end-2,2:end-1,2:end-1) )/2 ...
     + metrics.ET_X.*( rhop(2:end-1,3:end,2:end-1) - rhop(2:end-1,1:end-2,2:end-1) )/2 ...
     + metrics.ZT_X.*( rhop(2:end-1,2:end-1,3:end) - rhop(2:end-1,2:end-1,1:end-2) )/2 );
 
d_2 = 1/metrics.J.* ...
     ( metrics.XI_Y.*( rhop(3:end,2:end-1,2:end-1) - rhop(1:end-2,2:end-1,2:end-1) )/2 ...
     + metrics.ET_Y.*( rhop(2:end-1,3:end,2:end-1) - rhop(2:end-1,1:end-2,2:end-1) )/2 ...
     + metrics.ZT_Y.*( rhop(2:end-1,2:end-1,3:end) - rhop(2:end-1,2:end-1,1:end-2) )/2 );
 
d_3 = 1/metrics.J.* ...
     ( metrics.XI_Z.*( rhop(3:end,2:end-1,2:end-1) - rhop(1:end-2,2:end-1,2:end-1) )/2 ...
     + metrics.ET_Z.*( rhop(2:end-1,3:end,2:end-1) - rhop(2:end-1,1:end-2,2:end-1) )/2 ...
     + metrics.ZT_Z.*( rhop(2:end-1,2:end-1,3:end) - rhop(2:end-1,2:end-1,1:end-2) )/2 );

chi = 2*kappa*mean(d_1.^2 + d_2.^2 + d_3.^2,3);

end