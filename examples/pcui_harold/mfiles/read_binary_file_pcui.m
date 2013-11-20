function [x,y,z] =  read_binary_file_pcui(folder, filenameD, istep, params, isvector, isall)
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
ni = params.ni;
nj = params.nj;
nk = params.nk;

nni = ni/px+4;
nnj = nj/py+4;
nnk = nk/pz+4;

precision = 'float64';
N = nni*nnj*nnk;

if isvector
    ptrD = 8*(3*N+1);
else
    ptrD = 8*(N+1);
end

for indxz  = 1: pz
    for indxy = 1: py
        for indxx = 1: px               
            coord  = [ffD '.' num2str([2000+(indxx-1)*pz*py + (indxy-1)*pz + indxz-1])];
            fid = fopen(coord, 'r');
            fseek(fid, ptrD*(istep-1), 'bof');                    
            if isvector
                ddd = getdata(fid, 3*N, precision);
            else
                ddd = getdata(fid, N, precision);                   
            end
            fclose(fid);

            if isempty(ddd)
                break
            end

            for k = 1: nnk
                for j = 1: nnj
                    for i = 1: nni
                        if isvector
                            u(i, j, k) = ddd(nni*nnj*(k-1) + nni*(j-1) + i);
                            v(i, j, k) = ddd(N + nni*nnj*(k-1) + nni*(j-1) + i);
                            w(i, j, k) = ddd(2*N + nni*nnj*(k-1) +  nni*(j-1) + i);
                        else
                            s(i, j, k) = ddd(nni*nnj*(k-1) +  nni*(j-1) + i);
                        end
                    end
                end
            end

            if indxx == 1
                if isvector
                    if isall
                        U = u;
                        V = v;
                        W = w;
                    else
                        U = u(3:end-2, 3:end-2, 3:end-2);
                        V = v(3:end-2, 3:end-2, 3:end-2);
                        W = w(3:end-2, 3:end-2, 3:end-2);
                    end
                else    
                    if isall
                        S = s;
                    else                        
                        S = s(3:end-2, 3:end-2, 3:end-2);
                    end
                end
            else 
                if isvector
                    if isall
                        U = [U;u];
                        V = [V;v];
                        W = [W;w];
                    else
                        U = [U;u(3:end-2, 3:end-2, 3:end-2)];
                        V = [V;v(3:end-2, 3:end-2, 3:end-2)];
                        W = [W;w(3:end-2, 3:end-2, 3:end-2)];
                    end
                else
                    if isall
                        S = [S;s];
                    else
                        S = [S;s(3:end-2, 3:end-2, 3:end-2)];
                    end
                end
            end    

        end

        if isvector
            if isempty(ddd)
                break
            end
        end

        if indxy == 1
            if isvector
                UU = U;
                VV = V;
                WW = W;
            else
                SS = S;
            end
        else
            if isvector
                UU = cat(2, UU, U);
                VV = cat(2, VV, V);
                WW = cat(2, WW, W);
            else
                SS = cat(2, SS, S);
            end
        end
    end
    
    if indxz == 1
        if isvector
            UUU = UU;
            VVV = VV;
            WWW = WW;                
        else
            SSS = SS;
        end
    else
        if isvector
            UUU = cat(3, UUU, UU);
            VVV = cat(3, VVV, VV);
            WWW = cat(3, WWW, WW);
        else
            SSS = cat(3, SSS, SS);
        end
    end      
end

if isvector
    x = UUU;
    y = VVV;
    z = WWW;
else
    x = SSS;
end

