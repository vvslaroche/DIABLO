% script to load in full 3D snapshots and/or 2D slices for debugging 
% also for loading other variables e.g. S1, F1, etc
% specifically being used during duct.f timestepping

clear,clc%,close all

datapath = '/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-data/';
dirname = 'duct-knownmodes03-ai2';


fU1 = 1;
cU1 = ["9"];


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


%% plot U1,U2,U3,TH
% if fU1==1
%     for i=length(cU1)
%         a
%     end
% end

M = readPHYSICAL(NX,NYP,NZP,'U1_CHECK8',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'U1 CHECK8')

M = readPHYSICAL(NX,NYP,NZP,'U2_CHECK8',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'U2 CHECK8')

M = readPHYSICAL(NX,NYP,NZP,'U3_CHECK8',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'U3 CHECK8')

M = readPHYSICAL(NX,NYP,NZP,'TH_CHECK8',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'TH CHECK8')

%% R terms
M = readPHYSICAL(NX,NYP,NZP,'R1_CHECK7',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'R1 CHECK7')

M = readPHYSICAL(NX,NYP,NZP,'R2_CHECK7',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'R2 CHECK7')

M = readPHYSICAL(NX,NYP,NZP,'R3_CHECK7',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'R3 CHECK7')

M = readPHYSICAL(NX,NYP,NZP,'RTH_CHECK7',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'RTH CHECK7')

%% in fourier space
M = readFOURIER(NXP,NYP,NZ,'CU1_CHECK3',datapath,dirname);
plotFOURIER(M,NZ,'CU1 CHECK3')

M = readFOURIER(NXP,NYP,NZ,'CU2_CHECK3',datapath,dirname);
plotFOURIER(M,NZ,'CU2 CHECK3')

M = readFOURIER(NXP,NYP,NZ,'CU3_CHECK3',datapath,dirname);
plotFOURIER(M,NZ,'CU3 CHECK3')

M = readFOURIER(NXP,NYP,NZ,'CTH_CHECK3',datapath,dirname);
plotFOURIER(M,NZ,'CTH CHECK3')



%% for checking before and after REM_DIV (for comparison between separate runs)
M = readFOURIER(NXP,NYP,NZ,'CU2_CHECK2',datapath,dirname);
plotFOURIER(M,NZ,'CU2 CHECK2')
clim([-2.5e-3 2.5e-3])
M = readFOURIER(NXP,NYP,NZ,'CU3_CHECK2',datapath,dirname);
plotFOURIER(M,NZ,'CU3 CHECK2')
clim([-2.5e-3 2.5e-3])

M = readFOURIER(NXP,NYP,NZ,'CU2_CHECK3',datapath,dirname);
plotFOURIER(M,NZ,'CU2 CHECK3')
clim([-2.5e-3 2.5e-3])
M = readFOURIER(NXP,NYP,NZ,'CU3_CHECK3',datapath,dirname);
plotFOURIER(M,NZ,'CU3 CHECK3')
clim([-2.5e-3 2.5e-3])


%% for checking SAVE_STATS (applied BCs and FFT) after REM_DIV (for comparison between separate runs)
M = readFOURIER(NXP,NYP,NZ,'CU2_CHECK4',datapath,dirname);
plotFOURIER(M,NZ,'CU2 CHECK4')
clim([-2.5e-3 2.5e-3])
M = readFOURIER(NXP,NYP,NZ,'CU3_CHECK4',datapath,dirname);
plotFOURIER(M,NZ,'CU3 CHECK4')
clim([-2.5e-3 2.5e-3])

%% for checking CU at start of RK STEP
M = readFOURIER(NXP,NYP,NZ,'CU1_CHK1_RK1',datapath,dirname);
plotFOURIER(M,NZ,'CU1 CHECK1 RK1')
clim([-2.5e-3 2.5e-3])
M = readFOURIER(NXP,NYP,NZ,'CU2_CHK1_RK1',datapath,dirname);
plotFOURIER(M,NZ,'CU2 CHECK1 RK1')
clim([-2.5e-3 2.5e-3])
M = readFOURIER(NXP,NYP,NZ,'CU3_CHK1_RK1',datapath,dirname);
plotFOURIER(M,NZ,'CU3 CHECK1 RK1')
clim([-2.5e-3 2.5e-3])
M = readFOURIER(NXP,NYP,NZ,'CTH_CHK1_RK1',datapath,dirname);
plotFOURIER(M,NZ,'CTH CHECK1 RK1')


%% pressure
M = readFOURIER(NXP,NYP,NZ,'CP_CHK1_RK3',datapath,dirname);
plotFOURIER(M-max(M,[],'all'),NZ,'CP CHECK1 RK3')

temp = M(:,2:end,:,:) - M(:,1:end-1,:,:);
dMdz = zeros(size(M));
dMdz(:,2:end,:,:) = temp;
plotFOURIER(dMdz,NZ,'dCPdz')

temp = M(:,:,2:end,:) - M(:,:,1:end-1,:);
dMdy = zeros(size(M));
dMdy(:,:,2:end,:) = temp;
plotFOURIER(dMdy,NZ,'dCPdy')

% M = readFOURIER(NXP,NYP,NZ,'CR1_DUCT4',datapath,dirname);
% plotFOURIER(M,NZ,'CR1 DUCT4')
M = readFOURIER(NXP,NYP,NZ,'CR1_CHK9_RK2',datapath,dirname);
plotFOURIER(abs(M-max(M,[],'all')),NZ,'CR1 CHECK9 RK2')
clim([0 0.0005])

%% for checking U at end of RK STEP
M = readPHYSICAL(NX,NYP,NZP,'U1_CHK9_RK3',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'U1 CHECK9 RK3')
clim([-2.5e-3 2.5e-3])
M = readPHYSICAL(NX,NYP,NZP,'U2_CHK9_RK3',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'U2 CHECK9 RK3')
clim([-2.5e-3 2.5e-3])
M = readPHYSICAL(NX,NYP,NZP,'U3_CHK9_RK3',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'U3 CHECK9 RK3')
clim([-2.5e-3 2.5e-3])
M = readPHYSICAL(NX,NYP,NZP,'TH_CHK9_RK3',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'TH CHECK9 RK3')
% clim([-2.5e-3 2.5e-3])

%% for checking CU at end of RK STEP, before and after REM_DIV
M = readFOURIER(NXP,NYP,NZ,'CU1_CHK5_RK2',datapath,dirname);
plotFOURIER(M,NZ,'CU1 CHECK5 RK2')
clim([-2.5e-3 2.5e-3])
M = readFOURIER(NXP,NYP,NZ,'CU2_CHK5_RK2',datapath,dirname);
plotFOURIER(M,NZ,'CU2 CHECK5 RK2')
clim([-2.5e-3 2.5e-3])
M = readFOURIER(NXP,NYP,NZ,'CU3_CHK5_RK2',datapath,dirname);
plotFOURIER(M,NZ,'CU3 CHECK5 RK2')
clim([-2.5e-3 2.5e-3])
M = readFOURIER(NXP,NYP,NZ,'CTH_CHK5_RK2',datapath,dirname);
plotFOURIER(M,NZ,'CTH CHECK5 RK2')


M = readFOURIER(NXP,NYP,NZ,'CU1_CHK9_RK2',datapath,dirname);
plotFOURIER(M,NZ,'CU1 CHECK9 RK2')
clim([-2.5e-3 2.5e-3])
M = readFOURIER(NXP,NYP,NZ,'CU2_CHK9_RK2',datapath,dirname);
plotFOURIER(M,NZ,'CU2 CHECK9 RK2')
clim([-2.5e-3 2.5e-3])
M = readFOURIER(NXP,NYP,NZ,'CU3_CHK9_RK2',datapath,dirname);
plotFOURIER(M,NZ,'CU3 CHECK9 RK2')
clim([-2.5e-3 2.5e-3])
M = readFOURIER(NXP,NYP,NZ,'CTH_CHK9_RK2',datapath,dirname);
plotFOURIER(M,NZ,'CTH CHECK9 RK2')

%% checking R terms
M = readPHYSICAL(NX,NYP,NZP,'R1_CHK5_RK2',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'R1 CHECK5 RK2')
M = readPHYSICAL(NX,NYP,NZP,'R2_CHK5_RK2',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'R2 CHECK5 RK2')
M = readPHYSICAL(NX,NYP,NZP,'R3_CHK5_RK2',datapath,dirname);
plotPHYSICAL(M,NY,NYP,NZ,NZP,'R3 CHECK5 RK2')
% M = readPHYSICAL(NX,NYP,NZP,'RTH_CHK5_RK2',datapath,dirname);
% plotPHYSICAL(M,NY,NYP,NZ,NZP,'RTH CHECK5 RK2')

%% checking CR terms
M = readFOURIER(NXP,NYP,NZ,'CR1_CHK2_RK3',datapath,dirname);
plotFOURIER(M,NZ,'CR1 CHECK2 RK3')
M = readFOURIER(NXP,NYP,NZ,'CR2_CHK2_RK3',datapath,dirname);
plotFOURIER(M,NZ,'CR2 CHECK2 RK3')
M = readFOURIER(NXP,NYP,NZ,'CR3_CHK2_RK3',datapath,dirname);
plotFOURIER(M,NZ,'CR3 CHECK2 RK3')

%% checking CF terms
M = readFOURIER(NXP,NYP,NZ,'CF1_CHK3_RK1',datapath,dirname);
plotFOURIER(M,NZ,'CF1 CHECK3 RK1')
M = readFOURIER(NXP,NYP,NZ,'CF2_CHK3_RK1',datapath,dirname);
plotFOURIER(M,NZ,'CF2 CHECK3 RK1')
M = readFOURIER(NXP,NYP,NZ,'CF3_CHK3_RK1',datapath,dirname);
plotFOURIER(M,NZ,'CF3 CHECK3 RK1')


%% function definitions


function M = readPHYSICAL(NX,NYP,NZP,name,datapath,dirname)
M = nan(NX+2,NZP+2,NYP+2,4);
for i=1:4
    fname = [datapath dirname '/' name '_' char(i-1+48) '.dat'];
    temp = readlines(fname);
    data = zeros((NX+2)*(NZP+2)*(NYP+2),1);
    for j=1:length(data)
        if contains(temp(j),"*")
            data(j)=NaN;
        elseif contains(temp(j),"-Infinity")
            data(j)=-Inf;
        elseif contains(temp(j),"Infinity")
            data(j)=Inf;
%         elseif temp(j)==""
%             continue
        else
            data(j)=str2num(temp(j));
        end
    end
    M(:,:,:,i) = reshape(data,[NX+2 NZP+2 NYP+2]);
end
end


function plotPHYSICAL(M,NY,NYP,NZ,NZP,name)
temp1 = nan(2*(NYP+2)+1, 2*(NZP+2)+1);
temp1(1:NYP+2    ,1:NZP+2)     = squeeze(M(1,:,:,1))';
temp1(NYP+2+2:end,1:NZP+2)     = squeeze(M(1,:,:,2))';
temp1(1:NYP+2    ,NZP+2+2:end) = squeeze(M(1,:,:,3))';
temp1(NYP+2+2:end,NZP+2+2:end) = squeeze(M(1,:,:,4))';

temp2 = nan(NY, NZ);
temp2(1:NYP,1:NZP) = temp1(2:NYP+1,2:NZP+1);
temp2(NYP+1:end,1:NZP) = temp1(NYP+1+3+2:end-1,2:NZP+1);
temp2(1:NYP,NZP+1:end) = temp1(2:NYP+1,NZP+1+3+2:end-1);
temp2(NYP+1:end,NZP+1:end) = temp1(NYP+1+3+2:end-1,NZP+1+3+2:end-1);

figure
imagesc(temp1')
shading flat; colorbar; axis equal
xlabel('Y'), ylabel('Z')
title([name ' with processor gaps'])

figure
imagesc(temp2')
shading flat; colorbar; axis equal
xlabel('Y'), ylabel('Z')
title([name ' merged/trimmed'])
end



function M = readFOURIER(NXP,NYP,NZ,name,datapath,dirname)
M = nan(NXP+1,NZ+2,NYP+2,4);
for i=1:4
    fname = [datapath dirname '/' name '_' char(i-1+48) '.dat'];
    temp = readlines(fname);
    data = zeros((NXP+1)*(NZ+2)*(NYP+2),1);
    for j=1:length(data)
        if contains(temp(j),"*")
            data(j)=NaN;
        elseif contains(temp(j),"-Infinity")
            data(j)=-Inf;
        elseif contains(temp(j),"Infinity")
            data(j)=Inf;
%         elseif temp(j)==""
%             continue
        else
            data(j)=str2num(temp(j));
        end
    end
    M(:,:,:,i) = reshape(data,[NXP+1 NZ+2 NYP+2]);
end
end


function plotFOURIER(M,NZ,name)
figure
imagesc([squeeze(M(1,:,:,1)) nan(NZ+2,1) squeeze(M(1,:,:,2))])
shading flat; colorbar; axis equal
xlabel('Y'), ylabel('Z')
title([name ' with processor gaps'])
end