[uh,vh,wh] = read_binary_file_pcui(working_folder, 'output_UVW', istep, ...
                                      params, 1,1);
[xh,yh,zh] = read_binary_file_pcui(working_folder, 'xyz', 1, params,1,1);
                                  
nni = params.ni/params.px+4;
nnj = params.nj/params.py+4;
nnk = params.nk/params.pz+4;
exi = []; exj = []; exk = [];
if params.px == 1
    exi = 1:nni;
else
    for i=0:params.px-1
        if i==0
            exi = [exi,1:nni-2];
        elseif i==params.px-1
            exi = [exi,i*nni+3:(i+1)*nni];
        else
            exi = [exi,i*nni+3:(i+1)*nni-2];
        end
    end
end

if params.py == 1
    exj = 1:nnj;
else
    for j=0:params.py-1
        if j==0
            exj = [exj,1:nnj-2];
        elseif j==params.py-1
            exj = [exj,j*nnj+3:(j+1)*nnj];
        else
            exj = [exj,j*nnj+3:(j+1)*nnj-2];
        end
    end
end

if params.pz == 1
    exk = 1:nnk;
else
    for k=0:params.pz-1
        if k==0
            exk = [exk,1:nnk-2];
        elseif k==params.pz-1
            exk = [exk,k*nnk+3:(k+1)*nnk];
        else
            exk = [exk,k*nnk+3:(k+1)*nnk-2];
        end
    end
end


uh = uh(exi,exj,exk);
vh = vh(exi,exj,exk);
wh = wh(exi,exj,exk);
yh = yh(exi,exj,exk);

% imagesc(squeeze(uh(:,:,1))'); 
% colorbar; 
% colormap gray;
% 
figure;
plot(squeeze(uh(1,:,1)),squeeze(yh(1,:,1)),'ko-'); 
hold on; 
plot(squeeze(uh(2,:,1)),squeeze(yh(2,:,1)),'bo-'); 
plot(squeeze(uh(3,:,1)),squeeze(yh(3,:,1)),'mo-'); 
plot(squeeze(uh(4,:,1)),squeeze(yh(4,:,1)),'ro-');
load fastfit %load fast data
u_west = [ftanh(squeeze(yh(1,squeeze(yh(1,:,1))<jhoverlap,1))); ...
    flog(squeeze(yh(1,squeeze(yh(1,:,1))>=jhoverlap,1)))];
% plot(squeeze(uh(5,:,1)),squeeze(yh(5,:,1)),'co-');
% u_west = 1/6*squeeze(yh(1,:,1));
plot(u_west,squeeze(yh(1,:,1)),'k--','linewidth',2);
legend('n=-1','n=0','n=1','n=2','prescribed u_{west}','location','southeast');
xlabel('u [m/s]');
ylabel('z [m]');
% 
figure;
plot(squeeze(uh(end,:,1)),squeeze(yh(end,:,1)),'k<-'); 
hold on; 
plot(squeeze(uh(end-1,:,1)),squeeze(yh(end-1,:,1)),'b*'); 
plot(squeeze(uh(end-2,:,1)),squeeze(yh(end-2,:,1)),'mo'); 
plot(squeeze(uh(end-3,:,1)),squeeze(yh(end-3,:,1)),'co');
plot(squeeze(uh(end-4,:,1)),squeeze(yh(end-4,:,1)),'ro-');
legend('n=N+2','n=N+1','n=N','n=N-1','n=N-2','location','southeast');
xlabel('u [m/s]');
ylabel('z [m]');
% u_west = 1/6*squeeze(yh(1,:,1));
% % plot(u_west,squeeze(yh(1,:,1)),'g--');

figure;
plot(mean(uh(4,:,:),3),squeeze(yh(1,:,Ny/2)),'ko-','linewidth',2); 
hold on; 
plot(mean(uh(Nx/4,:,:),3),squeeze(yh(Nx/4,:,Ny/2)),'bo-','linewidth',2); 
plot(mean(uh(Nx/2,:,:),3),squeeze(yh(Nx/2,:,Ny/2)),'mo-','linewidth',2); 
plot(mean(uh(3*Nx/4,:,:),3),squeeze(yh(3*Nx/4,:,Ny/2)),'co-','linewidth',2);
plot(mean(uh(Nx-2,:,:),3),squeeze(yh(Nx,:,Ny/2)),'ro-','linewidth',2);
legend('x/L=0','x/L=0.25','x/L=0.5','x/L=0.75','x/L=1','location','northwest')
ylim([0 0.3])
xlabel('u [m s^{-1}]','fontsize',14)
ylabel('z [m]','fontsize',14)
% print('-f4','-r500','-dpng','development')
u_west = 1/6*squeeze(yh(1,:,1));
% plot(u_west,squeeze(yh(1,:,1)),'g--');

figure;
plot(squeeze(wh(1,:,3)),squeeze(yh(1,:,1)),'ko-'); 
hold on; 
plot(squeeze(wh(2,:,3)),squeeze(yh(2,:,1)),'bo-'); 
plot(squeeze(wh(3,:,3)),squeeze(yh(3,:,1)),'mo-'); 
plot(squeeze(wh(4,:,3)),squeeze(yh(4,:,1)),'ro-');
load fastfit %load fast data
u_west = [ftanh(squeeze(yh(1,squeeze(yh(1,:,1))<jhoverlap,1))); ...
    flog(squeeze(yh(1,squeeze(yh(1,:,1))>=jhoverlap,1)))];
% plot(squeeze(uh(5,:,1)),squeeze(yh(5,:,1)),'co-');
% u_west = 1/6*squeeze(yh(1,:,1));
% plot(u_west,squeeze(yh(1,:,1)),'k--','linewidth',2);
legend('n=-1','n=0','n=1','n=2','location','southeast');
xlabel('w [m/s]');
ylabel('z [m]');