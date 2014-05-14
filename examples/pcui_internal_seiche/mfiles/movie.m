clc; clear all; close all;
load 512x512_midslice.mat

%%
figure;
hold on;
for n=1:1
    display(n);
    cla;
    pcolor(squeeze(x),squeeze(z),squeeze(w(:,:,20)));
    shading flat;
    axis image;
    colorbar;
    drawnow;
%     pause;
end