function xpart =  read_binary_particles_pcui(folder, filenameD, istep, params)
%
% Filename    : read_binary_particles_pcui.m

% Author      : Bobby Arthur

% Description : Reads in particle  data output by PCUI. Expects unformatted 
%               Fortran binary files. PArticle data is only output on
%               processor 0. Data is stored in matlab vectors.

ffD = fullfile(folder, filenameD);

precision = 'float64';

ni = params.ni;
nj = params.nj;
nk = params.nk;
% npart = (ni+2)*(nj+2)*(nk+2);
npart = ni*nj;
N = npart*3;
ptrD = 8*N;
              
fid = fopen(ffD, 'r');
fseek(fid, ptrD*(istep-1), 'bof');
ddd = getdata(fid, N, precision);
fclose(fid);

% keyboard;
if isempty(ddd)
    return
end

%extract particle values
xpart = zeros(npart,3);
xpart(:,1) = ddd(1:npart);
xpart(:,2) = ddd(npart+1:2*npart);
xpart(:,3) = ddd(2*npart+1:3*npart);

end