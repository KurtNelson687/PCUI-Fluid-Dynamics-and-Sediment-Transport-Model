function [ddd,Ep,Eb,dEbdt,phi_d] =  read_binary_energy_pcui(folder, filenameD, dt)
%
% Filename    : read_binary_energy_pcui.m

% Author      : Bobby Arthur

% Description : Reads in energy data output by PCUI. Expects unformatted 
%               Fortran binary files. Energy data is only output on
%               processor 0. Energy values are store in matlab vectors.

ffD = fullfile(folder, filenameD);

precision = 'float64';
              
coord  = [ffD '.500'];
fid = fopen(coord, 'r');
ddd = getdata(fid, inf, precision);
fclose(fid);

if isempty(ddd)
    return
end

nout = 3;
Ep = ddd(1:nout+1:end);
Eb = ddd(2:nout+1:end);
dEbdt = (Eb(3:end)-Eb(1:end-2))/2/dt;
phi_d = ddd(3:nout+1:end);

end

