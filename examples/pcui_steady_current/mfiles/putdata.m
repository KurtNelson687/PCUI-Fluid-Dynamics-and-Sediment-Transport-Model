%function data = getdata(fid, arraysize, precision)
function data = putdata(fid, array, precision)

  %
  % SGI is fopen(fname,'rb,'ieee-be');
  % Linux is fopen(fname,'rb','native');
  %

  tic;
  if(fid>0),
%     if(precision=='float32'),
%       space = 4*(arraysize + 1);                    
%     elseif(precision=='float64')
%       space = 8*(arraysize + 1);                    
%     end
    headerflag = 'float32'; % This needs to be changed to float64 to read in C binary data
    footerflag = headerflag;

    a = fwrite(fid,1,headerflag);%header
    data = fwrite(fid,array,precision);%data
    b = fwrite(fid,1, footerflag);%footer

  else
    fprintf('Error in function: getdata...undefined fid\n');
    data = zeros(array,1);
  end
  
%  fprintf('Read data at a rate of %.2f Mb/sec\n',space/1e6/toc);