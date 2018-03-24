function [a] =  read_binary_file_pcui_profile(folder, filenameD, istep, params, isall,numGhost)
%
% Filename    : read_binary_file_pcui.m

% Author      : Goncalo Gil

% Description : Reads in data output by PCUI. Expects unformatted Fortran
%               binary files. It reads each file created by a particular
%               processor and assembles the data into a single Matlab
%               array.
%
%             - If isvector is equal to 0 it reads in three-dimensional
%               arrays (i.e. density field)
%
%             - If isvector is equal to 1 it reads in four-dimensional
%               arrays (i.e. grid and velocity field field)
%
%             - The istep variable indicates the time step at which the
%               information is extracted. If the grid (x,y,z) is requested,
%               then set istep = 1.
%
%             - If isall = 1 then it extracts all the ghost points for each
%               submatrix. If isall = 0 then it only extracts the values of
%               the array at the physical grid points.
%
% Suggested improvements: 1) Vectorize operations in loop
%                         2) Keep at least one ghost point of from exterior
%                            processors in order to interpolate value to
%                            boundary. (PCUI does not have a grid point on
%                            the boundary)
%
% Author      : Goncalo Gil, Stanford University
% email       : gilg@stanford.edu
%

ffD = fullfile(folder, filenameD);

px = params.px;
py = params.py;
pz = params.pz;
nj = params.nj;

nnj = nj/py+2*numGhost;

precision = 'float64';
N = nnj;

ptrD = 8*(N+1);

for indxy = 1: py
    coord  = [ffD '.' num2str(100+ (indxy-1))];
    fid = fopen(coord, 'r');
    fseek(fid, ptrD*(istep-1), 'bof');
    ddd = getdata(fid, N, precision);
    fclose(fid);
    
    if indxy == 1
        a = ddd;
    else
        a = [a; ddd];
    end
end
end
