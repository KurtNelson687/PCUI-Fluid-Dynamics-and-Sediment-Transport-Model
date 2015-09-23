function metrics_2d = calculate_binary_2d_metrics(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates metric quantities in 2 dimensions (x,y)-(horizontal, vertical)
% using binary output from PCUI. x,y are the grid points (cell centers) 
% including 1 halo cell.
%
% Bobby Arthur
% March 2015
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_xi = 0.5*( x(3:end,2:end-1) - x(1:end-2,2:end-1) );
Y_xi = 0.5*( y(3:end,2:end-1) - y(1:end-2,2:end-1) );
X_et = 0.5*( x(2:end-1,3:end) - x(2:end-1,1:end-2) );
Y_et = 0.5*( y(2:end-1,3:end) - y(2:end-1,1:end-2) );

metrics_2d.J = X_xi.*Y_et - X_et.*Y_xi;

metrics_2d.XI_X =  Y_et;
metrics_2d.ET_X = -Y_xi;
metrics_2d.XI_Y = -X_et;
metrics_2d.ET_Y =  X_xi;

end