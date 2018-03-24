function write_binary_file_pcui_profile(folder, filenameX, params, A)
precision = 'float64';

fn = fullfile(folder,filenameX);
ff  = [fn '.1'];
fid = fopen(ff, 'w');
putdata(fid,A,precision);
fclose(fid);
end