%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script creates a movie showing the time evolution of
% velocity magnitude. It loads the data file "velPlanes.mat"
% which is created by post processing PCUI output files with 
% "ExtractHPCData.mat".
%
% Author      : Kurt Nelson, Stanford University
% email       : knelson3@stanford.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
working_folder = '/Users/kurtnelson/Desktop/HPC_transfer/videoForOliver/';
load([working_folder 'velPlanes.mat'])
moviename = 'TurbSpinUp.gif';
frameTime =0.2;
skipFrames = 1;


[m,n,p] = size(velMag);
fig1 = figure(1);
hold

count = 1;
frames = 1:skipFrames:800;
for numFrame = frames
    
    pcolor(x_xyplane,y_xyplane,velMag(:,:,numFrame));
    xlab = xlabel('$x$');
    ylab = ylabel('$z$');
    set(xlab,'interpreter','Latex','FontSize',12)
    set(ylab,'interpreter','Latex','FontSize',12)
    
    %axis([0 max(max(x_plotvar1)) 0 max(max(x_plotvar2))])
    %set(gca,'YTick',[min(min(x_plotvar2)) .05 0.1],'YTickLabel',{'0', '0.5', '0.1'})
    axis image
    shading flat
    %colorbar
    %caxis([0 0.2])
    
    
    frame = getframe(gcf);
    if count == 1
        [im,map] = rgb2ind(frame.cdata,256,'nodither');
        im(1,1,1,length(frames)) = 0;
    else
        im(:,:,1,count) = rgb2ind(frame.cdata,map,'nodither');
    end
    clf
    count=count+1;
end

imwrite(im,map,[working_folder '/' moviename],'DelayTime',frameTime,'LoopCount',inf);

