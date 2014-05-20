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
yh = yh(exi,exj,exk);

imagesc(squeeze(uh(:,:,1))'); 
colorbar; 
colormap gray;

figure;
plot(squeeze(uh(1,:,1)),squeeze(yh(1,:,1)),'ko-'); 
hold on; 
plot(squeeze(uh(2,:,1)),squeeze(yh(2,:,1)),'bo-'); 
plot(squeeze(uh(3,:,1)),squeeze(yh(3,:,1)),'mo-'); 
plot(squeeze(uh(4,:,1)),squeeze(yh(4,:,1)),'ro-');
% plot(squeeze(uh(5,:,1)),squeeze(yh(5,:,1)),'co-');
u_west = 1/6*squeeze(yh(1,:,1));
plot(u_west,squeeze(yh(1,:,1)),'k--');
legend('n=-1','n=0','n=1','n=2','prescribed u_{west}','location','southeast');
xlabel('u [m/s]');
ylabel('z [m]');

figure;
plot(squeeze(uh(end,:,1)),squeeze(yh(end,:,1)),'ko-'); 
hold on; 
plot(squeeze(uh(end-1,:,1)),squeeze(yh(end-1,:,1)),'bo-'); 
plot(squeeze(uh(end-2,:,1)),squeeze(yh(end-2,:,1)),'mo'); 
plot(squeeze(uh(end-3,:,1)),squeeze(yh(end-3,:,1)),'co-');
plot(squeeze(uh(end-4,:,1)),squeeze(yh(end-4,:,1)),'r--');
u_west = 1/6*squeeze(yh(1,:,1));
plot(u_west,squeeze(yh(1,:,1)),'g--');

figure;
plot(squeeze(uh(1,:,1)),squeeze(yh(1,:,1)),'ko-'); 
hold on; 
plot(squeeze(uh(32,:,1)),squeeze(yh(32,:,1)),'bo-'); 
plot(squeeze(uh(64,:,1)),squeeze(yh(64,:,1)),'mo'); 
plot(squeeze(uh(96,:,1)),squeeze(yh(96,:,1)),'co-');
plot(squeeze(uh(128,:,1)),squeeze(yh(128,:,1)),'r--');
u_west = 1/6*squeeze(yh(1,:,1));
plot(u_west,squeeze(yh(1,:,1)),'g--');
