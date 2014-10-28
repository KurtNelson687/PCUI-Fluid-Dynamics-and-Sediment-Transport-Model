function [x,y,z] =  read_binary_particle_grid_pcui(folder, filenameD, istep, params)
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

ni = params.ni+2;
nj = params.nj+2;
nk = params.nk+2;

precision = 'float64';
N = ni*nj*nk;
ptrD = 8*(3*N+1);
           
fid = fopen(ffD, 'r');
fseek(fid, ptrD*(istep-1), 'bof');                    
ddd = getdata(fid, 3*N, precision);
fclose(fid);

for k = 1:nk
    for j = 1:nj
        for i = 1:ni
                x(i, j, k) = ddd(ni*nj*(k-1) + ni*(j-1) + i);
                y(i, j, k) = ddd(N + ni*nj*(k-1) + ni*(j-1) + i);
                z(i, j, k) = ddd(2*N + ni*nj*(k-1) +  ni*(j-1) + i);
        end
    end
end

end

