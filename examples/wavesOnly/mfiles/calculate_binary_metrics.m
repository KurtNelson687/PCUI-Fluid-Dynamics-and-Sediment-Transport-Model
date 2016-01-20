function metrics = calculate_binary_metrics(x,y,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates metric quantities using binary output from cns-koltakov. x,y,z 
% are the grid points (cell centers) including 1 halo cell.
%
% Bobby Arthur
% May 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_xi = 0.5*( x(3:end,2:end-1,2:end-1) - x(1:end-2,2:end-1,2:end-1) );
Y_xi = 0.5*( y(3:end,2:end-1,2:end-1) - y(1:end-2,2:end-1,2:end-1) );
Z_xi = 0.5*( z(3:end,2:end-1,2:end-1) - z(1:end-2,2:end-1,2:end-1) );

X_et = 0.5*( x(2:end-1,3:end,2:end-1) - x(2:end-1,1:end-2,2:end-1) );
Y_et = 0.5*( y(2:end-1,3:end,2:end-1) - y(2:end-1,1:end-2,2:end-1) );
Z_et = 0.5*( z(2:end-1,3:end,2:end-1) - z(2:end-1,1:end-2,2:end-1) );

X_zt = 0.5*( x(2:end-1,2:end-1,3:end) - x(2:end-1,2:end-1,1:end-2) );
Y_zt = 0.5*( y(2:end-1,2:end-1,3:end) - y(2:end-1,2:end-1,1:end-2) );
Z_zt = 0.5*( z(2:end-1,2:end-1,3:end) - z(2:end-1,2:end-1,1:end-2) );

metrics.J = ( X_xi.*Y_et.*Z_zt + X_et.*Y_zt.*Z_xi + X_zt.*Y_xi.*Z_et ...
    - X_zt.*Y_et.*Z_xi - X_et.*Y_xi.*Z_zt - X_xi.*Y_zt.*Z_et );

metrics.XI_X = ( Y_et.*Z_zt - Y_zt.*Z_et );
metrics.ET_X = ( Y_zt.*Z_xi - Y_xi.*Z_zt );
metrics.ZT_X = ( Y_xi.*Z_et - Y_et.*Z_xi );
metrics.XI_Y = ( Z_et.*X_zt - Z_zt.*X_et );
metrics.ET_Y = ( Z_zt.*X_xi - Z_xi.*X_zt );
metrics.ZT_Y = ( Z_xi.*X_et - Z_et.*X_xi );
metrics.XI_Z = ( X_et.*Y_zt - X_zt.*Y_et );
metrics.ET_Z = ( X_zt.*Y_xi - X_xi.*Y_zt );
metrics.ZT_Z = ( X_xi.*Y_et - X_et.*Y_xi );

end