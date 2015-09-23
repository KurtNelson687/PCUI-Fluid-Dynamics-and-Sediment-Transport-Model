function [P] = calculate_binary_production(up,vp,wp, ...
                                                 ubar,vbar,wbar,metrics_2d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates production of TKE using binary output from PCUI.
% up,vp,wp are fluctuating velocities with one halo cell.
% ubar,vbar,wbar are mean (laterally-averaged) velocities with 
% one halo cell.
%
% Bobby Arthur
% March 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate Reynolds Stresses
r_11 = -mean(up.*up,3);
r_22 = -mean(vp.*vp,3);
% r_33 = -mean(wp.*wp,3);
r_12 = -mean(up.*vp,3);
% r_13 = -mean(up.*wp,3);
% r_23 = -mean(vp.*wp,3);

%Calculate components of mean strain rate S = dU_i/dx_j
S_11 = 1/metrics_2d.J.* ...
     ( metrics_2d.XI_X.*( ubar(3:end,2:end-1,2:end-1) - ubar(1:end-2,2:end-1,2:end-1) )/2 ...
     + metrics_2d.ET_X.*( ubar(2:end-1,3:end,2:end-1) - ubar(2:end-1,1:end-2,2:end-1) )/2);
 
S_12 = 1/metrics_2d.J.* ...
     ( metrics_2d.XI_Y.*( ubar(3:end,2:end-1,2:end-1) - ubar(1:end-2,2:end-1,2:end-1) )/2 ...
     + metrics_2d.ET_Y.*( ubar(2:end-1,3:end,2:end-1) - ubar(2:end-1,1:end-2,2:end-1) )/2);
 
S_21 = 1/metrics_2d.J.* ...
     ( metrics_2d.XI_X.*( vbar(3:end,2:end-1,2:end-1) - vbar(1:end-2,2:end-1,2:end-1) )/2 ...
     + metrics_2d.ET_X.*( vbar(2:end-1,3:end,2:end-1) - vbar(2:end-1,1:end-2,2:end-1) )/2);
 
S_22 = 1/metrics_2d.J.* ...
     ( metrics_2d.XI_Y.*( vbar(3:end,2:end-1,2:end-1) - vbar(1:end-2,2:end-1,2:end-1) )/2 ...
     + metrics_2d.ET_Y.*( vbar(2:end-1,3:end,2:end-1) - vbar(2:end-1,1:end-2,2:end-1) )/2);

% S_31 = 1/metrics_2d.J.* ...
%      ( metrics_2d.XI_X.*( wbar(3:end,2:end-1,2:end-1) - wbar(1:end-2,2:end-1,2:end-1) )/2 ...
%      + metrics_2d.ET_X.*( wbar(2:end-1,3:end,2:end-1) - wbar(2:end-1,1:end-2,2:end-1) )/2);
%  
% S_32 = 1/metrics_2d.J.* ...
%      ( metrics_2d.XI_Y.*( wbar(3:end,2:end-1,2:end-1) - wbar(1:end-2,2:end-1,2:end-1) )/2 ...
%      + metrics_2d.ET_Y.*( wbar(2:end-1,3:end,2:end-1) - wbar(2:end-1,1:end-2,2:end-1) )/2);

%Calculate components of production
P_11 = -r_11(2:end-1,2:end-1).*S_11;
 
P_12 = -0.5*r_12(2:end-1,2:end-1).*(S_12 + S_21);
 
% P_13 = -0.5*r_13(2:end-1,2:end-1).*(S_13 + S_31);
 
P_22 = -r_22(2:end-1,2:end-1).*S_22;
 
% P_23 = -0.5*r_23(2:end-1,2:end-1).*(S_23 + S_32);

% P_33 = -r_33(2:end-1,2:end-1).*S_33;

% P = P_11 + 2*P_12 + 2*P_13 + P_22 + 2*P_23 + P_33;

P = P_11 + 2*P_12 + P_22;

end
