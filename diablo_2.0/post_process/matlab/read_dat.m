clear, close all

NX = 128;
NY = 33;
NZ = 17;
NZP = 9;

NPROCS = 4;
NPROCY = 2;
NPROCZ = NPROCS/NPROCY;


NXP = NX/2/NPROCZ;

% savepath = ['/Users/vincent/Documents/_School/Berkeley/Research/' ...
%     'DIABLO-data/duct-chan-knownmodes/'];
savepath = ['/Users/vincent/Documents/_School/Berkeley/Research/' ...
    'DIABLO-data/duct-knownmodes/'];


fname = [savepath 'U1_0.dat'];
f = fopen(fname);
data = fscanf(f,'%f');
fclose(f);
U1 = reshape(data,[NX+2 NZP+2 NY+2]);

fname = [savepath 'CU1_0.dat'];
f = fopen(fname);
data = fscanf(f,'%f');
fclose(f);
CU1 = reshape(data,[NXP+1 NZ+2 NY+2]);

fname = [savepath 'U1_new_0.dat'];
f = fopen(fname);
data = fscanf(f,'%f');
fclose(f);
U1_new = reshape(data,[NX+2 NZP+2 NY+2]);

fname = [savepath 'TMP1A_0.dat'];
f = fopen(fname);
data = fscanf(f,'%f');
fclose(f);
TMP1A = reshape(data,[NX/2 NZP+2 NY+2]);

fname = [savepath 'TMP2_0.dat'];
f = fopen(fname);
data = fscanf(f,'%f');
fclose(f);
TMP2 = reshape(data,[NXP (NZP+2)*NPROCZ]);

fname = [savepath 'TMP2_new_0.dat'];
f = fopen(fname);
data = fscanf(f,'%f');
fclose(f);
TMP2_new = reshape(data,[NXP (NZP+2)*NPROCZ]);

%% Plotting
cl = 0.01;

figure
imagesc(U1(:,:,5)'), colorbar
clim([-1 1])
axis equal, axis tight
set(gcf, 'Position',[100 400 1300 300])

figure
imagesc(U1_new(:,:,5)'), colorbar
clim([-1 1])
axis equal, axis tight
set(gcf, 'Position',[100 400 1300 300])
% return

% figure
% imagesc(TMP1A(:,:,5)'), colorbar
% clim([-cl cl])
% axis equal, axis tight
% set(gcf, 'Position',[100 400 1300 300])

figure
imagesc(TMP2(:,:)'), colorbar
clim([-cl cl])
axis equal, axis tight
set(gcf, 'Position',[100 400 1300 300])

figure
imagesc(TMP2_new(:,:)'), colorbar
clim([-cl cl])
axis equal, axis tight
set(gcf, 'Position',[100 400 1300 300])

% figure
% temp=TMP1A(1:32,:,5)-TMP2(:,1:11);
% imagesc(temp'), colorbar
% clim([-cl cl])
% axis equal, axis tight
% set(gcf, 'Position',[100 400 1300 300])


figure
% CU1(:,7:11,5)=CU1(:,7:11,5)*NZ;
imagesc(CU1(:,:,5)'), colorbar
clim([-cl cl])
axis equal, axis tight
set(gcf, 'Position',[100 400 1300 300])


% figure
% imagesc(squeeze(CU1(:,5,:))'), colorbar
% clim([-cl cl])
% axis equal, axis tight
% set(gcf, 'Position',[100 400 1300 300])




%% junk
% figure
% imagesc(squeeze(TMP1A(:,5,:))'), colorbar
% clim([-cl cl])
% axis equal, axis tight
% set(gcf, 'Position',[100 400 1300 300])

% data = fread(f,[128+2 8+2]);
% [32+1 8+2]
% data = fread(f,[128+2 8+2],'real*8');
% data = fread(f,[128+2 8+2],'float');
% data = readmatrix(savepath);
% data = fread(f,[32+1 8+2],'real*8');

% pcolor(data'), shading flat, colorbar
% imagesc(data'./max(max(data))), colorbar

% CU1 = reshape(data,[NXP+1 NZ+2]);
