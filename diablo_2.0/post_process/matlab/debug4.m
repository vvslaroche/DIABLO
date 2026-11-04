% this script reads in MAT?_C_YZ files to explore matrix solve

clc,clear

datapath = '/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-data/';
dirname = 'duct-knownmodes03';

NX = 16;
NY = 65;
NYP = 33;
NZ = 49;
NZP = 25;

NPROCZ = 2;
NXP = NX/(2*NPROCZ);

filename = [datapath dirname '/grid.h5'];
gy = h5read(filename,'/grids/y');
gz = h5read(filename,'/grids/z');
gyf = (gy(1:end-1)+gy(2:end))/2;
gzf = (gz(1:end-1)+gz(2:end))/2;
dy = diff(gyf);
dz = diff(gzf);
dyf = diff(gy);
dzf = diff(gz);


%% 
% output of matrices and vector
MATL_C_YZ = nan(NYP+2,NZ+2,4);
MATD_C_YZ = nan(NYP+2,NZ+2,4);
MATU_C_YZ = nan(NYP+2,NZ+2,4);
VEC_C_YZ = nan(NYP+2,NZ+2,4);


for i=1:4
%     fname = [datapath dirname '/MATL_C_YZ_' char(i-1+48) '.dat'];
%     f = fopen(fname);
%     data = fscanf(f,'%f');
%     fclose(f);
%     MATL_C_YZ(:,:,i) = reshape(data,[NYP+2 NZ+2]);
% 
%     fname = [datapath dirname '/MATD_C_YZ_' char(i-1+48) '.dat'];
%     f = fopen(fname);
%     data = fscanf(f,'%f');
%     fclose(f);
%     MATD_C_YZ(:,:,i) = reshape(data,[NYP+2 NZ+2]);
% 
%     fname = [datapath dirname '/MATU_C_YZ_' char(i-1+48) '.dat'];
%     f = fopen(fname);
%     data = fscanf(f,'%f');
%     fclose(f);
%     MATU_C_YZ(:,:,i) = reshape(data,[NYP+2 NZ+2]);

    fname = [datapath dirname '/VEC_C_YZ_CHECK9_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    VEC_C_YZ(:,:,i) = reshape(data,[NYP+2 NZ+2]);
end

figure
imagesc([squeeze(VEC_C_YZ(:,:,1))' nan(NZ+2,1) squeeze(VEC_C_YZ(:,:,2))'])
colorbar
axis equal
xlabel('Y'), ylabel('Z')
% title('VEC\_C\_YZ before TDMA')
% title('VEC new BC Y check 7')
title('VEC\_C\_YZ check 9')

%% assemble into big matrix - Z SWEEP
MAT = zeros((NYP+2)*(NZ+2),(NYP+2)*(NZ+2));
VEC = zeros((NYP+2)*(NZ+2),1);

for j=1:NYP+2
    for k=1:NZ+2
        MAT((j-1)*(NZ+2)+k,(j-1)*(NZ+2)+k) = MATD_C_YZ(j,k,1);
        if (j*k)<((NYP+2)*(NZ+2))
            MAT((j-1)*(NZ+2)+k,(j-1)*(NZ+2)+k+1) = MATU_C_YZ(j,k,1);
        end
        if (j*k)>(1)
            MAT((j-1)*(NZ+2)+k,(j-1)*(NZ+2)+k-1) = MATL_C_YZ(j,k,1);
        end
        VEC((j-1)*(NZ+2)+k) = VEC_C_YZ(j,k,1);
    end
end

SOL = MAT\VEC;

% alternate version that solves subset only (J=1,NY and K=1,NZ)
MAT = zeros((NYP)*(NZ),(NYP)*(NZ));
VEC = zeros((NYP)*(NZ),1);
for j=1:NYP
    for k=1:NZ
        MAT((j-1)*(NZ)+k,(j-1)*(NZ)+k) = MATD_C_YZ(j+1,k+1,1);
        if (j*k)<((NYP)*(NZ))
            MAT((j-1)*(NZ)+k,(j-1)*(NZ)+k+1) = MATU_C_YZ(j+1,k+1,1);
        end
        if (j*k)>(1)
            MAT((j-1)*(NZ)+k,(j-1)*(NZ)+k-1) = MATL_C_YZ(j+1,k+1,1);
        end
        VEC((j-1)*(NZ)+k) = VEC_C_YZ(j+1,k+1,1);
    end
end
SOL = MAT\VEC;

SOL_C_YZ = nan(NYP+2,NZ+2,4);
for j=1:NYP
    for k=1:NZ
        SOL_C_YZ(j+1,k+1,1) = SOL((j-1)*(NZ)+k);
    end
end


% SOL_C_YZ = nan(NYP+2,NZ+2,4);
% for j=1:NYP+2
%     for k=1:NZ+2
%         SOL_C_YZ(j,k,1) = SOL((j-1)*(NZ+2)+k);
%     end
% end

figure
imagesc([squeeze(SOL_C_YZ(:,:,1))' nan(NZ+2,1) squeeze(SOL_C_YZ(:,:,2))'])
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('SOL\_C\_YZ')


%% load in CU2 and CU3 from just before REM_DIV and calculate the first VEC
CU2_CHECK = nan(NXP+1,NZ+2,NYP+2,4);
CU3_CHECK = nan(NXP+1,NZ+2,NYP+2,4);
for i=1:4
    fname = [datapath dirname '/CU2_CHECK2_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    CU2_CHECK(:,:,:,i) = reshape(data,[NXP+1 NZ+2 NYP+2]);

    fname = [datapath dirname '/CU3_CHECK2_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    CU3_CHECK(:,:,:,i) = reshape(data,[NXP+1 NZ+2 NYP+2]);
end

figure
imagesc([squeeze(CU2_CHECK(1,:,:,1)) nan(NZ+2,1) squeeze(CU2_CHECK(1,:,:,2))])
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('CU2 with gaps')

figure
imagesc([squeeze(CU3_CHECK(1,:,:,1)) nan(NZ+2,1) squeeze(CU3_CHECK(1,:,:,2))])
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('CU3 with gaps')


VEC = zeros(NZ+2,NYP+2,2);
for j=2:NYP+1
    for k=2:NZ+1
%         VEC(k,j,1) = (CU2_CHECK(1,k,j+1,1)-CU2_CHECK(1,k,j,1))/dyfp(j) ...
%             +        (CU3_CHECK(1,k+1,j,1)-CU3_CHECK(1,k,j,1))/dzf(k);
% %         VEC(k,j,1) = (CU3_CHECK(1,k+1,j,1)-CU3_CHECK(1,k,j,1))/dzf(k);
% %         VEC(k,j,2) = (CU2_CHECK(1,k,j+1,2)-CU2_CHECK(1,k,j,2))/dyfp(j) ...
% %             +        (CU3_CHECK(1,k+1,j,2)-CU3_CHECK(1,k,j,2))/dzf(k);
    end
end

figure
imagesc([squeeze(VEC(:,:,1)) nan(NZ+2,1) squeeze(VEC(:,:,2))])
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('VEC calculated, with gaps')
clim([-0.01 0.01])

% function i=jk2i(j,k,NZ)
% i=(j-1)*NZ+k;
% end