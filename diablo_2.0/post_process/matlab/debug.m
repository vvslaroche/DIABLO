% script to load in full 3D snapshots and/or 2D slices for debugging 

clear,clc%,close all

datapath = '/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-data/';
dirname = 'duct-knownmodes03';
% dirname = 'dwl01';



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

LX = 30;
gx = linspace(0,LX,NX);
% LZ = 5;
% gzf = linspace(0,LZ,size(U,3));






%% reading full 3D output from a standard DIABLO 'out' file
% filename = [datapath dirname '/end.h5'];
filename = [datapath dirname '/output/out00000.h5'];
U = h5read(filename,'/Timestep/U');
V = h5read(filename,'/Timestep/V');
W = h5read(filename,'/Timestep/W');





figure
% imagesc(gyf,gzf,squeeze(V(1,:,:))')
imagesc(squeeze(V(1,:,:))')
% shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('V out')
% clim([-1 1])

figure
% imagesc(gyf,gzf,squeeze(W(1,:,:))')
imagesc(squeeze(W(1,:,:))')
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('W out')
% clim([-1 1])



%% reading 2D slice from movie file
filename = [datapath dirname '/movie.h5'];
U_zy = h5read(filename,'/u_zy/0001');
V_zy = h5read(filename,'/v_zy/0001');
W_zy = h5read(filename,'/w_zy/0001');
TH_zy = h5read(filename,'/th01_zy/0002');

figure
% imagesc(gyf,gzf,V_zy)
imagesc(V_zy)
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('V movie')
% clim([-1 1])

figure
% imagesc(gyf,gzf,W_zy)
imagesc(W_zy)
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('W movie')
% clim([-1 1])



return

%% reading 3D output via dat text files - PHYSICAL
U2_CHECK = nan(NX+2,NZP+2,NYP+2,4);
for i=1:4
    fname = [datapath dirname '/U3_CHECK9_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    U2_CHECK(:,:,:,i) = reshape(data,[NX+2 NZP+2 NYP+2]);
end


temp = nan(2*(NYP+2)+1, 2*(NZP+2)+1);
temp(1:NYP+2,1:NZP+2)         = squeeze(U2_CHECK(1,:,:,1))';
temp(NYP+2+2:end,1:NZP+2)     = squeeze(U2_CHECK(1,:,:,2))';
temp(1:NYP+2,NZP+2+2:end)     = squeeze(U2_CHECK(1,:,:,3))';
temp(NYP+2+2:end,NZP+2+2:end) = squeeze(U2_CHECK(1,:,:,4))';

figure
imagesc(temp')
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('V with gaps')

temp2 = nan(NY, NZ);
temp2(1:NYP,1:NZP) = temp(2:NYP+1,2:NZP+1);
temp2(NYP+1:end,1:NZP) = temp(NYP+1+3+2:end-1,2:NZP+1);
temp2(1:NYP,NZP+1:end) = temp(2:NYP+1,NZP+1+3+2:end-1);
temp2(NYP+1:end,NZP+1:end) = temp(NYP+1+3+2:end-1,NZP+1+3+2:end-1);

figure
imagesc(temp2')
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('V merged/trimmed')


%% reading 3D output via dat text files - FOURIER
CU2_CHECK = nan(NXP+1,NZ+2,NYP+2,4);
for i=1:4
    fname = [datapath dirname '/CU3_CHECK3_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    CU2_CHECK(:,:,:,i) = reshape(data,[NXP+1 NZ+2 NYP+2]);
end


% temp = nan(2*(NYP+2)+1, 2*(NZP+2)+1);
% temp(1:NYP+2,1:NZP+2)         = squeeze(U2_CHECK(1,:,:,1))';
% temp(NYP+2+2:end,1:NZP+2)     = squeeze(U2_CHECK(1,:,:,2))';
% temp(1:NYP+2,NZP+2+2:end)     = squeeze(U2_CHECK(1,:,:,3))';
% temp(NYP+2+2:end,NZP+2+2:end) = squeeze(U2_CHECK(1,:,:,4))';

figure
% imagesc(temp')
imagesc([squeeze(CU2_CHECK(1,:,:,1)) nan(NZ+2,1) squeeze(CU2_CHECK(1,:,:,2))])
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('CU3 CHECK 3 with gaps')
% clim([-3 3])

% temp2 = nan(NY, NZ);
% temp2(1:NYP,1:NZP) = temp(2:NYP+1,2:NZP+1);
% temp2(NYP+1:end,1:NZP) = temp(NYP+1+3+2:end-1,2:NZP+1);
% temp2(1:NYP,NZP+1:end) = temp(2:NYP+1,NZP+1+3+2:end-1);
% temp2(NYP+1:end,NZP+1:end) = temp(NYP+1+3+2:end-1,NZP+1+3+2:end-1);
% 
% figure
% imagesc(temp2')
% shading flat
% colorbar
% axis equal
% xlabel('Y'), ylabel('Z')
% title('CU2 merged/trimmed')


%% compare before and after REM_DIV
CU_CHECK2 = nan(NXP+1,NZ+2,NYP+2,4);
CU_CHECK3 = nan(NXP+1,NZ+2,NYP+2,4);
for i=1:4
    fname = [datapath dirname '/CU2_CHECK2_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    CU_CHECK2(:,:,:,i) = reshape(data,[NXP+1 NZ+2 NYP+2]);

    fname = [datapath dirname '/CU2_CHECK3_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    CU_CHECK3(:,:,:,i) = reshape(data,[NXP+1 NZ+2 NYP+2]);
end
temp = CU_CHECK3-CU_CHECK2;

figure
imagesc([squeeze(temp(1,:,:,1)) nan(NZ+2,1) squeeze(temp(1,:,:,2))])
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('CU2 (CHECK 3 - CHECK 2) with gaps')



for i=1:4
    fname = [datapath dirname '/CU3_CHECK2_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    CU_CHECK2(:,:,:,i) = reshape(data,[NXP+1 NZ+2 NYP+2]);

    fname = [datapath dirname '/CU3_CHECK3_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    CU_CHECK3(:,:,:,i) = reshape(data,[NXP+1 NZ+2 NYP+2]);
end
temp = CU_CHECK3-CU_CHECK2;

figure
imagesc([squeeze(temp(1,:,:,1)) nan(NZ+2,1) squeeze(temp(1,:,:,2))])
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('CU3 (CHECK 3 - CHECK 2) with gaps')


%% reading files from the middle of FFT routines - TMP1A
TMP1A = nan(NX/2,NZP+2,NYP+2,4);
for i=1:4
    fname = [datapath dirname '/TMP1A_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    TMP1A(:,:,:,i) = reshape(data,[NX/2 NZP+2 NYP+2]);
end

figure
% imagesc(temp')
imagesc([ [squeeze(TMP1A(1,:,:,1)) nan(NZP+2,1) squeeze(TMP1A(1,:,:,2))]; ...
           nan(1,2*(NYP+2)+1); ...
          [squeeze(TMP1A(1,:,:,3)) nan(NZP+2,1) squeeze(TMP1A(1,:,:,4))] ])
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('TMP1A with gaps')
% clim([-0.00001 .00001])



%% reading files from the middle of FFT routines - TMP2
TMP2 = nan(NXP,(NZP+2)*NPROCZ,NYP+2,4);
for i=1:4
    fname = [datapath dirname '/TMP2_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    TMP2(:,:,:,i) = reshape(data,[NXP (NZP+2)*NPROCZ NYP+2]);
end

% figure
% % imagesc(temp')
% % imagesc([ [squeeze(TMP2(:,:,1,1)) nan((NZP+2)*NPROCZ,1) squeeze(TMP2(:,:,1,2))]; ...
% %            nan(1,2*(NYP+2)+1); ...
% %           [squeeze(TMP2(:,:,1,3)) nan((NZP+2)*NPROCZ,1) squeeze(TMP2(:,:,1,4))] ])
% imagesc([ [squeeze(TMP2(:,:,1,1))  squeeze(TMP2(:,:,1,2))]; ...
% %            nan(1,2*(NYP+2)+1); ...
%           [squeeze(TMP2(:,:,1,3))  squeeze(TMP2(:,:,1,4))] ])
% shading flat
% colorbar
% axis equal
% xlabel('X'), ylabel('Z')
% title('TMP1A with gaps')
% 
% figure
% imagesc(squeeze(TMP2(:,:,1,1))')
% axis equal
% xlabel('X'), ylabel('Z')
 

figure
imagesc([squeeze(TMP2(1,:,:,1)) nan((NZP+2)*NPROCZ,1) squeeze(TMP2(1,:,:,2))])
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('TMP2 with gaps')

%%
% output of matrices and vector
MATL_C_YZ = nan(NYP+2,NZ+2,4);
MATD_C_YZ = nan(NYP+2,NZ+2,4);
MATU_C_YZ = nan(NYP+2,NZ+2,4);
VEC_C_YZ = nan(NYP+2,NZ+2,4);

% fname = [datapath dirname '/MATD_C_YZ_0.dat'];
% f = fopen(fname);
% data = fscanf(f,'%f');
% fclose(f);
% MATD_C_YZ_0 = reshape(data,[NX+2 NZP+2 NYP+2]);

for i=1:4
    fname = [datapath dirname '/MATL_C_YZ_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    MATL_C_YZ(:,:,i) = reshape(data,[NYP+2 NZ+2]);

    fname = [datapath dirname '/MATD_C_YZ_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    MATD_C_YZ(:,:,i) = reshape(data,[NYP+2 NZ+2]);

    fname = [datapath dirname '/MATU_C_YZ_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    MATU_C_YZ(:,:,i) = reshape(data,[NYP+2 NZ+2]);

    fname = [datapath dirname '/VEC_C_YZ_' char(i-1+48) '.dat'];
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
title('VEC\_C\_YZ before TDMA')




%% from dat files - CU2 and CR1 after TDMA
CU2 = nan(NXP+1,NZ+2,NYP+2,4);
CR1 = nan(NXP+1,NZ+2,NYP+2,4);
for i=1:4
%     fname = [datapath dirname '/CU2_' char(i-1+48) '.dat'];
%     f = fopen(fname);
%     data = fscanf(f,'%f');
%     fclose(f);
%     CU2(:,:,:,i) = reshape(data,[NXP+1 NZ+2 NYP+2]);

    fname = [datapath dirname '/CR1_CHECK4_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    CR1(:,:,:,i) = reshape(data,[NXP+1 NZ+2 NYP+2]);
end


% figure
% imagesc([squeeze(CU2(1,:,:,1)) nan(NZ+2,1) squeeze(CU2(1,:,:,2))])
% colorbar
% axis equal
% xlabel('Y'), ylabel('Z')
% title('CU2 before adding CR1')

figure
imagesc([squeeze(CR1(1,:,:,1)) nan(NZ+2,1) squeeze(CR1(1,:,:,2))])
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('CR1 check 4')



% % output from start of save_stats_duct - full 3D U2
% fname = [datapath dirname '/U2_insavestats0.dat'];
% f = fopen(fname);
% data = fscanf(f,'%f');
% fclose(f);
% V_YZ_dat0 = reshape(data,[NX+2 NZP+2 NYP+2]);
% 
% fname = [datapath dirname '/U2_insavestats1.dat'];
% f = fopen(fname);
% data = fscanf(f,'%f');
% fclose(f);
% V_YZ_dat1 = reshape(data,[NX+2 NZP+2 NYP+2]);
% 
% fname = [datapath dirname '/U2_insavestats2.dat'];
% f = fopen(fname);
% data = fscanf(f,'%f');
% fclose(f);
% V_YZ_dat2 = reshape(data,[NX+2 NZP+2 NYP+2]);
% 
% fname = [datapath dirname '/U2_insavestats3.dat'];
% f = fopen(fname);
% data = fscanf(f,'%f');
% fclose(f);
% V_YZ_dat3 = reshape(data,[NX+2 NZP+2 NYP+2]);
% 
% 
% temp = nan(2*(NYP+2)+1, 2*(NZP+2)+1);
% temp(1:NYP+2,1:NZP+2)         = squeeze(V_YZ_dat0(1,:,:))';
% temp(NYP+2+2:end,1:NZP+2)     = squeeze(V_YZ_dat1(1,:,:))';
% temp(1:NYP+2,NZP+2+2:end)     = squeeze(V_YZ_dat2(1,:,:))';
% temp(NYP+2+2:end,NZP+2+2:end) = squeeze(V_YZ_dat3(1,:,:))';
% 
% figure
% imagesc(temp')
% shading flat
% colorbar
% axis equal
% xlabel('Y'), ylabel('Z')
% title('V dat SSD')
% 
% V_YZ_dat = zeros(NY, NZ);
% V_YZ_dat(1:NYP,1:NZP) = temp(2:NYP+1,2:NZP+1);
% V_YZ_dat(NYP+1:end,1:NZP) = temp(NYP+1+3+2:end-1,2:NZP+1);
% V_YZ_dat(1:NYP,NZP+1:end) = temp(2:NYP+1,NZP+1+3+2:end-1);
% V_YZ_dat(NYP+1:end,NZP+1:end) = temp(NYP+1+3+2:end-1,NZP+1+3+2:end-1);
% 
% figure
% imagesc(gyf,gzf,V_YZ_dat')
% shading flat
% colorbar
% axis equal
% xlabel('Y'), ylabel('Z')
% title('V dat SSD')




% fname = [datapath dirname '/V_YZ_insavestats0.dat'];
% f = fopen(fname);
% data = fscanf(f,'%f');
% fclose(f);
% V_YZ_hdf0 = reshape(data,[NZP NYP]);
% 
% fname = [datapath dirname '/V_YZ_insavestats1.dat'];
% f = fopen(fname);
% data = fscanf(f,'%f');
% fclose(f);
% V_YZ_hdf1 = reshape(data,[NZP NYP]);
% 
% fname = [datapath dirname '/V_YZ_insavestats2.dat'];
% f = fopen(fname);
% data = fscanf(f,'%f');
% fclose(f);
% V_YZ_hdf2 = reshape(data,[NZP NYP]);
% 
% fname = [datapath dirname '/V_YZ_insavestats3.dat'];
% f = fopen(fname);
% data = fscanf(f,'%f');
% fclose(f);
% V_YZ_hdf3 = reshape(data,[NZP NYP]);
% 
% temp = nan(NY+1, NZ+1);
% temp(1:NYP,1:NZP)         = V_YZ_hdf0';
% temp(NYP+2:end,1:NZP)     = V_YZ_hdf1(:,2:end)';
% temp(1:NYP,NZP+2:end)     = V_YZ_hdf2(2:end,:)';
% temp(NYP+2:end,NZP+2:end) = V_YZ_hdf3(2:end,2:end)';
% 
% figure
% imagesc(temp')
% shading flat
% colorbar
% axis equal
% xlabel('Y'), ylabel('Z')
% title('V dat in hdf5')





% figure
% % imagesc(gyf,gzf,W_zy)
% imagesc(squeeze(V_YZ_dat(1,:,:)))
% shading flat
% colorbar
% axis equal
% % xlabel('Y'), ylabel('Z')
% title('V dat')



%% check divergence
dudx = (U(2:end,:,:) - U(1:end-1,:,:)) ./ reshape((gx(2:end)-gx(1:end-1)),[length(gx)-1,1,1]);
dvdy = (V(:,2:end,:) - V(:,1:end-1,:)) ./ reshape((gyf(2:end)-gyf(1:end-1)),[1,length(gyf)-1,1]);
dwdz = (W(:,:,2:end) - W(:,:,1:end-1)) ./ reshape((gzf(2:end)-gzf(1:end-1)),[1,1,length(gzf)-1]);

div = sum(dudx,'all') + sum(dvdy,'all') + sum(dwdz,'all');





% %% function
% 
% function [xp,yp,Ap] = prepcolor(x,y,A)
%     xp = zerosx;
%     xp(end+1) = 0;
% 
% end