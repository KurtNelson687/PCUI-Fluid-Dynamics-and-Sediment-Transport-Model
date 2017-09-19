function du = fftdiff(u,L,m)

     sizeu = size(u);
     if(sizeu(1)==1)
       u=u';
     end

     N = length(u);
     u = reshape(u,N,1);
     k = 2*pi*[-N/2:N/2-1]'/L;
     uhat = fftshift(fft(u));
     duhat = (1i*k).^m.*uhat;
     du = real(ifft(fftshift(duhat)));

     if(sizeu(1)==1)
       du=du';
     end
end

