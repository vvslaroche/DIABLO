% script to look at CHANNEL version
% check initial condition
% explore REM_DIV subroutine

clc,clear

datapath = '/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-data/';
dirname = 'duct-chan-knownmodes03';

% NX = 16;
% NY = 65;
% NYP = 33;
% NZ = 48;
% NZP = 24;

NX = 48;
NY = 65;
NYP = 33;
NZ = 16;
NZP = 8;

NPROCZ = 2;
NXP = NX/(2*NPROCZ);

filename = [datapath dirname '/grid.h5'];
gy = h5read(filename,'/grids/y');
gyf = (gy(1:end-1)+gy(2:end))/2;
% gzf = (gz(1:end-1)+gz(2:end))/2;
% dy = diff(gyf);
% dz = diff(gzf);
% dyf = diff(gy);
% dzf = diff(gz);


%% reading 3D output via dat text files - PHYSICAL
U_CHECK = nan(NX+2,NZP+2,NYP+2,4);
for i=1:4
    fname = [datapath dirname '/U1_CHECK1_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    U_CHECK(:,:,:,i) = reshape(data,[NX+2 NZP+2 NYP+2]);
end

figure
imagesc([ [squeeze(U_CHECK(1,:,:,1)) nan(NZP+2,1) squeeze(U_CHECK(1,:,:,2))]; ...
           nan(1,2*(NYP+2)+1); ...
          [squeeze(U_CHECK(1,:,:,3)) nan(NZP+2,1) squeeze(U_CHECK(1,:,:,4))] ])
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('U1 CHECK 1 with gaps')


%% reading 3D output via dat text files - PHYSICAL XY
U_CHECK = nan(NX+2,NZP+2,NYP+2,4);
for i=1:4
    fname = [datapath dirname '/U1_CHECK2_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    U_CHECK(:,:,:,i) = reshape(data,[NX+2 NZP+2 NYP+2]);
end

figure
imagesc([squeeze(U_CHECK(:,1,:,1)) nan(NX+2,1) squeeze(U_CHECK(:,1,:,2))])
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('X')
title('U1 CHECK 2 with gaps')


%% reading 3D output via dat text files - FOURIER
CU_CHECK = nan(NXP+1,NZ+2,NYP+2,4);
for i=1:4
    fname = [datapath dirname '/CU2_CHECK2_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    CU_CHECK(:,:,:,i) = reshape(data,[NXP+1 NZ+2 NYP+2]);
end

figure
imagesc([squeeze(CU_CHECK(1,:,:,1)) nan(NZ+2,1) squeeze(CU_CHECK(1,:,:,2))])
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('CU2 CHECK 2 with gaps')
clim([-0.5 0.5])


%% reading full 3D output from a standard DIABLO 'out' file
% filename = [datapath dirname '/end.h5'];
filename = [datapath dirname '/output/out00002.h5'];
U = h5read(filename,'/Timestep/U');
V = h5read(filename,'/Timestep/V');
W = h5read(filename,'/Timestep/W');


figure
imagesc(squeeze(V(1,:,:))')
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('V out')
% clim([-1 1])

figure
imagesc(squeeze(W(1,:,:))')
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('W out')
% clim([-1 1])


%% read VEC
VEC_C = nan(NXP+1,NYP+2,4);


for i=1:4
    fname = [datapath dirname '/VEC_C_CHECK3_' char(i-1+48) '.dat'];
    f = fopen(fname);
    data = fscanf(f,'%f');
    fclose(f);
    VEC_C(:,:,i) = reshape(data,[NXP+1 NYP+2]);
end

figure
imagesc([ [squeeze(VEC_C(:,:,1)) nan(NXP+1,1) squeeze(VEC_C(:,:,2))]; ...
           nan(1,2*(NYP+2)+1); ...
          [squeeze(VEC_C(:,:,3)) nan(NXP+1,1) squeeze(VEC_C(:,:,4))] ])
% imagesc([ [squeeze(U_CHECK(1,:,:,1)) nan(NZP+2,1) squeeze(U_CHECK(1,:,:,2))]; ...
%            nan(1,2*(NYP+2)+1); ...
%           [squeeze(U_CHECK(1,:,:,3)) nan(NZP+2,1) squeeze(U_CHECK(1,:,:,4))] ])
colorbar
axis equal
xlabel('Y'), ylabel('X')
title('VEC\_C check 3')