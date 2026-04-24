% this is a matlab implementation of a parallel LU decomposition solver
% goal is to practice constructing the variables and assembling into full
% matrix 
% variables have to be in complex space - so that Z is not parallelized

clc,clear,close all


NY = 10;
NZ = 8;
LY = 50;
LZ = 30;
GY = linspace(0,LY,NY);
% GY = linspace(1,1.5,NY).*GY;
GZ = linspace(0,LZ,NZ);
% GZ = linspace(1,1.5,NZ).*GZ;
DY = (diff(GY(1:end-1))+diff(GY(2:end)))/2;
DY = [DY(1) DY DY(end)];
DZ = (diff(GZ(1:end-1))+diff(GZ(2:end)))/2;
DZ = [DZ(1) DZ DZ(end)];


MATD_C_YZ   = zeros(NY,NZ);
MATL_C_YZ   = zeros(NY,NZ);
MATLL_C_YZ  = zeros(NY,NZ);
MATU_C_YZ   = zeros(NY,NZ);
MATUU_C_YZ  = zeros(NY,NZ);
VEC_C_YZ    = zeros(NY,NZ);
MAT_C_YZ    = zeros(NY*NZ,2*NZ+1);
% could potentially format as 
% MAT_C_YZ_alt    = zeros(NY,NZ,2*NZ+1);

M = zeros(NY*NZ,NY*NZ);


for j=1:NY
    for k=1:NZ
        disp(['element j,k=' num2str(j) ',' num2str(k) ' is number ' num2str((j-1)*NZ+k)])

        
        MATD_C_YZ(j,k) = 1/(DY(j))^2 + 1/(DY(j))^2  +  1/(DZ(k))^2 + 1/(DZ(k))^2;
        if j>1
            MATLL_C_YZ(j,k) = -1/(DY(j-1))^2;
        end
        if j<NY
            MATUU_C_YZ(j,k) = -1/(DY(j+1))^2;
        end
        if k>1
            MATL_C_YZ(j,k) = -1/(DZ(k-1))^2;
        end
        if k<NZ
            MATU_C_YZ(j,k) = -1/(DZ(k+1))^2;
        end


        MAT_C_YZ((j-1)*NZ+k,NZ+1) = 1/(DY(j))^2 + 1/(DY(j))^2  +  1/(DZ(k))^2 + 1/(DZ(k))^2;
        if j>1
            MAT_C_YZ((j-1)*NZ+k,1) = -1/(DY(j-1))^2;
        end
        if j<NY
            MAT_C_YZ((j-1)*NZ+k,2*NZ+1) = -1/(DY(j+1))^2;
        end
        if k>1
            MAT_C_YZ((j-1)*NZ+k,NZ+1-1) = -1/(DZ(k-1))^2;
        end
        if k<NZ
            MAT_C_YZ((j-1)*NZ+k,NZ+1+1) = -1/(DZ(k+1))^2;
        end


        % comparison full matrix construction
        M((j-1)*NZ+k,(j-1)*NZ+k) = M((j-1)*NZ+k,(j-1)*NZ+k)+ 1/(DY(j))^2 + 1/(DY(j))^2  +  1/(DZ(k))^2 + 1/(DZ(k))^2;
        if j>1
            M((j-1)*NZ+k,(j-1-1)*NZ+k) = M((j-1)*NZ+k,(j-1-1)*NZ+k) -1/(DY(j-1))^2;
%             M((j-1)*NZ+k,(j-1)*NZ+k-1) = M((j-1)*NZ+k,(j-1-1)*NZ+k) -1/(DY(j-1))^2;
        end
        if j<NY
            M((j-1)*NZ+k,(j-1+1)*NZ+k) = M((j-1)*NZ+k,(j-1+1)*NZ+k) -1/(DY(j+1))^2;
%             M((j-1)*NZ+k,(j-1)*NZ+k+1) = M((j-1)*NZ+k,(j-1+1)*NZ+k) -1/(DY(j+1))^2;
        end
        if k>1
            M((j-1)*NZ+k,(j-1)*NZ+k-1) = M((j-1)*NZ+k,(j-1)*NZ+k-1) -1/(DZ(k-1))^2;
        end
        if k<NZ
            M((j-1)*NZ+k,(j-1)*NZ+k+1) = M((j-1)*NZ+k,(j-1)*NZ+k+1) -1/(DZ(k+1))^2;
        end
    end
end


return
[L,U] = lu(M);
% x =

% forward sweep
% calculate L and U and store in MAT_C_YZ
% also solve for y of the system
% Ax=b splits into: Ly=b and Ux=y
for i1=1:NY*NZ
    for i2=1:2*NZ+1
        % calculate L and put into left side of MAT
        if i2<NZ+1 && i2>NZ+1-i1
            MAT_C_YZ(i1,i2) = MAT_C_YZ(i1,i2);
            for i3=1:i2-1
                MAT_C_YZ(i1,i2) = MAT_C_YZ(i1,i2) - U(i3,i2)*L(i1,i3);
            end
            MAT_C_YZ(i1,i2) = MAT_C_YZ(i1,i2)/MAT_C_YZ(i1-(NZ+1-i2),NZ+1);
        end

        % calculate U and put into right side of MAT
        if i2>=NZ+1
            MAT_C_YZ(i1,i2) = MAT_C_YZ(i1,i2);
            for i3=1:i1-1
                MAT_C_YZ(i1,i2) = MAT_C_YZ(i1,i2) - MAT_C_YZ(i3,i2)*L(i1,i3);
            end
        end
    end



%     L(i1,i1) = 1;
% 
%     lowerbound = 1;
%     if i1>NZ+1
%         lowerbound = i1-NZ;
%     end
%     upperbound = NY*NZ;
%     if i1<NY*NZ-NZ
%         upperbound = i1+NZ;
%     end
% 
%     for i2=lowerbound:upperbound
%         if i1>i2
%             L(i1,i2) = A(i1,i2)/U(i2,i2);
%             for i3=1:i2-1
%                 L(i1,i2) = L(i1,i2) - U(i3,i2)*L(i1,i3)/U(i2,i2);
%                 ops = ops+1;
%             end
%         end
% 
%         if i2>=i1
%             U(i1,i2) = A(i1,i2);
%             for i3=1:i1-1
%                 U(i1,i2) = U(i1,i2) - U(i3,i2)*L(i1,i3);
%                 ops = ops+1;
%             end
%         end
%     end


    if i1==25; break; end
end