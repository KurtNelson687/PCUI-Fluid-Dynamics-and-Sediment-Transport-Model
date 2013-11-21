function write_binary_file_pcui(folder, filenameX, params, A)

fn = fullfile(folder,filenameX);
px = params.px;
py = params.py;
pz = params.pz;
nni = params.ni/params.px+4;
nnj = params.nj/params.py+4;
nnk = params.nk/params.pz+4;

precision = 'float64';

for indxz  = 1: pz
    for indxy = 1: py
        for indxx = 1: px  
            ff  = [fn '.' num2str(700+(indxx-1)*pz*py + (indxy-1)*pz + indxz-1)];
            fid = fopen(ff, 'w');
            putdata(fid,A((indxx-1)*nni+1:indxx*nni, ...
                          (indxy-1)*nnj+1:indxy*nnj, ...
                          (indxz-1)*nnk+1:indxz*nnk,:),precision);
            fclose(fid);
        end
    end
end