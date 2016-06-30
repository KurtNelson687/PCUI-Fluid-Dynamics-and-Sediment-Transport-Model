%
% Filename    : assemble_array.m
% Author      : Goncalo Gil
% Description : Assembles a Matlab array (3-D or 4-D) into a PCUI ready set
%               of binary files depending on the number of processors
%               requested. Essentially it takes a singular Matlab array
%               with no ghostpoints and creates a series of subarrays with 
%               ghostpoints as required by PCUI.
%                                       
% Author      : Goncalo Gil, Stanford University
% email       : gilg@stanford.edu
%

function A = assemble_array(B,params)

px = params.px;
py = params.py;
pz = params.pz;
nni_A = params.ni/params.px+4;  
nnj_A = params.nj/params.py+4;
nnk_A = params.nk/params.pz+4;
nni_B = nni_A-4;
nnj_B = nnj_A-4;
nnk_B = nnk_A-4;
xb_s = 0;
yb_s = 0;
zb_s = 0;
xb_e = 0;
yb_e = 0;
zb_e = 0;
nd = length(B(1,1,1,:));
periodicx = 1;
periodicy = 1;
periodicz = 1;
A = zeros(nni_A*px,nnj_A,nnk_A,nd);

for indxz  = 1: pz
    for indxy = 1: py
        for indxx = 1: px
           
            % Interior (i)
            % X direction -> B = Back; F = Front
            ix_A = (indxx-1) * nni_A + 3:indxx * nni_A - 2;
            ix_B = (indxx-1) * nni_B + 1:indxx * nni_B;                         
            B_A = (indxx-1) * nni_A + 1:(indxx-1) * nni_A + 2;
            F_A = indxx * nni_A - 1:indxx * nni_A;
            
            % Deal with special cases
            if indxx == 1 && indxx == px                
                B_B = [2 1];                
                F_B = [indxx * nni_B indxx * nni_B - 1];
            elseif indxx == 1
                B_B = [2 1];                
                F_B = [indxx * nni_B + 1 indxx * nni_B + 2];
            elseif indxx == px
                B_B = [(indxx-1) * nni_B - 1 (indxx-1) * nni_B];
                F_B = [indxx * nni_B indxx * nni_B - 1];
            else
                B_B = (indxx-1) * nni_B - 1:(indxx-1) * nni_B;
                F_B = indxx * nni_B + 1:indxx * nni_B + 2;
            end
            
            % Y direction -> L = Left; R = Right 
            iy_A = (indxy-1) * nnj_A + 3:indxy * nnj_A - 2;
            iy_B = (indxy-1) * nnj_B + 1:indxy * nnj_B;                       
            L_A = (indxy-1) * nnj_A + 1:(indxy-1) * nnj_A + 2;
            R_A = indxy * nnj_A - 1:indxy * nnj_A;
            
            % Deal with special cases
            if indxy == 1 && indxy == py                
                L_B = [2 1];                
                R_B = [indxy * nnj_B indxy * nnj_B - 1];
            elseif indxy == 1
                L_B = [2 1];                
                R_B = [indxy * nnj_B + 1 indxy * nnj_B + 2];
            elseif indxy == py
                L_B = [(indxy-1) * nnj_B - 1 (indxy-1) * nnj_B];
                R_B = [indxy * nnj_B indxy * nnj_B - 1];
            else
                L_B = (indxy-1) * nnj_B - 1:(indxy-1) * nnj_B;
                R_B = indxy * nnj_B + 1:indxy * nnj_B + 2;              
            end           
            
            % Z direction -> D = Down; U = Up            
            iz_A = (indxz-1) * nnk_A + 3:indxz * nnk_A - 2;
            iz_B = (indxz-1) * nnk_B + 1:indxz * nnk_B;                                        
            D_A = (indxz-1) * nnk_A + 1:(indxz-1) * nnk_A + 2;
            U_A = indxz * nnk_A - 1:indxz * nnk_A;
            
         % Deal with special cases
            if indxz == 1 && indxz == pz                
                D_B = [2 1];                
                U_B = [indxz * nnk_B indxz * nnk_B - 1];
            elseif indxz == 1
                D_B = [2 1];
                U_B = [indxz * nnk_B + 1 indxz * nnk_B + 2];
            elseif indxz == pz
                D_B = [(indxz-1) * nnk_B - 1 (indxz-1) * nnk_B];
                U_B = [indxz * nnk_B indxz * nnk_B - 1];           
            else
                D_B = (indxz-1) * nnk_B - 1:(indxz-1) * nnk_B;
                U_B = indxz * nnk_B + 1:indxz * nnk_B + 2;              
            end           
            
            % Kernel (1)
            A(ix_A,iy_A,iz_A,:) = B(ix_B,iy_B,iz_B,:);
            
            % Corners (8)
            A(B_A,L_A,D_A,:) = B(B_B,L_B,D_B,:);
            A(F_A,L_A,D_A,:) = B(F_B,L_B,D_B,:);
            A(B_A,R_A,D_A,:) = B(B_B,R_B,D_B,:);
            A(F_A,R_A,D_A,:) = B(F_B,R_B,D_B,:);
            A(B_A,L_A,U_A,:) = B(B_B,L_B,U_B,:);
            A(F_A,L_A,U_A,:) = B(F_B,L_B,U_B,:);
            A(B_A,R_A,U_A,:) = B(B_B,R_B,U_B,:);
            A(F_A,R_A,U_A,:) = B(F_B,R_B,U_B,:);
            
            % Columns (12)
            A(ix_A,L_A,D_A,:) = B(ix_B,L_B,D_B,:);
            A(ix_A,L_A,U_A,:) = B(ix_B,L_B,U_B,:);
            A(ix_A,R_A,D_A,:) = B(ix_B,R_B,D_B,:);
            A(ix_A,R_A,U_A,:) = B(ix_B,R_B,U_B,:);
            A(B_A,iy_A,U_A,:) = B(B_B,iy_B,U_B,:);
            A(B_A,iy_A,D_A,:) = B(B_B,iy_B,D_B,:);
            A(F_A,iy_A,U_A,:) = B(F_B,iy_B,U_B,:);
            A(F_A,iy_A,D_A,:) = B(F_B,iy_B,D_B,:);
            A(B_A,L_A,iz_A,:) = B(B_B,L_B,iz_B,:);
            A(F_A,L_A,iz_A,:) = B(F_B,L_B,iz_B,:);
            A(B_A,R_A,iz_A,:) = B(B_B,R_B,iz_B,:);
            A(F_A,R_A,iz_A,:) = B(F_B,R_B,iz_B,:);
            
            % Sides (6)
            
            A(ix_A,iy_A,U_A,:) = B(ix_B,iy_B,U_B,:);
            A(ix_A,iy_A,D_A,:) = B(ix_B,iy_B,D_B,:);
            A(B_A,iy_A,iz_A,:) = B(B_B,iy_B,iz_B,:);
            A(F_A,iy_A,iz_A,:) = B(F_B,iy_B,iz_B,:);
            A(ix_A,L_A,iz_A,:) = B(ix_B,L_B,iz_B,:);
            A(ix_A,R_A,iz_A,:) = B(ix_B,R_B,iz_B,:);
            
            % If sub-array is on the boundary:
            if indxx == 1
                if periodicx
                    A(1,:,:,1) = A(end-3,:,:,1);
                    A(2,:,:,1) = A(end-2,:,:,1);
                else
                    A(1,:,:,1) = 2*xb_s-A(4,:,:,1);
                    A(2,:,:,1) = 2*xb_s-A(3,:,:,1);
                end
            end
            
            if indxy == 1 && nd > 1
                if periodicy
                    A(:,1,:,2) = A(:,end-3,:,2);
                    A(:,2,:,2) = A(:,end-2,:,2);
                else
                    A(:,1,:,2) = 2*yb_s-A(:,4,:,2);
                    A(:,2,:,2) = 2*yb_s-A(:,3,:,2);
                end
            end
            
            if indxz == 1 && nd > 2
                if periodicz
                    A(:,:,1,3) = A(:,:,end-3,3);
                    A(:,:,2,3) = A(:,:,end-2,3);
                else
                    A(:,:,1,3) = 2*zb_s-A(:,:,4,3);
                    A(:,:,2,3) = 2*zb_s-A(:,:,3,3);
                end
            end
            
            if indxx == px
                if periodicx
                    A(end-1,:,:,1) = A(1,:,:,1);
                    A(end,:,:,1) = A(2,:,:,1);
                else
                    A(end-1,:,:,1) = 2*xb_e-A(end-2,:,:,1);
                    A(end,:,:,1)   = 2*xb_e-A(end-3,:,:,1);
                end
            end
            
            if indxy == py && nd > 1
                if periodicy
                    A(:,end-1,:,2) = A(:,1,:,2);
                    A(:,end,:,2) = A(:,2,:,2);
                else
                    A(:,end-1,:,2) = 2*yb_e-A(:,end-2,:,2);
                    A(:,end,:,2)   = 2*yb_e-A(:,end-3,:,2);
                end
            end
            
            if indxz == pz && nd > 2
                if periodicz
                    A(:,:,end-1,3) = A(:,:,1,3);
                    A(:,:,end,3) = A(:,:,2,3);
                else
                    A(:,:,end-1,3) = 2*zb_e-A(:,:,end-2,3);
                    A(:,:,end,3)   = 2*zb_e-A(:,:,end-3,3);
                end
            end

        end
    end
end