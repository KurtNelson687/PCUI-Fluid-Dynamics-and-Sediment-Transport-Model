function [a] =  read_binary_file_pcui_value(folder, filenameD, istep, numpoints)
%
% Filename    : read_binary_file_pcui_value.m

% Author      : Modified verision of Goncalo Gil's read_binary_file_pcui - modified by Kurt Nelson

% Description : Reads in single point data from PCUI output. Expects unformatted Fortran
%               binary files. 

%             - The istep variable indicates the time step at which the
%               information is extracted. 
%              
%              - numpoints is the number of data points to be extracted at
%              each timestep.

ffD = fullfile(folder, filenameD);

nnj = numpoints;

precision = 'float64';
N = nnj;

ptrD = 8*(N+1);

for indxy = 1
    coord  = [ffD '.' num2str(100+ (indxy-1))];
    fid = fopen(coord, 'r');
    fseek(fid, ptrD*(istep-1), 'bof');
    ddd = getdata(fid, N, precision);
    fclose(fid);
    a = ddd;
end
end
