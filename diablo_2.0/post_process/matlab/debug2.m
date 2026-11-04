% for creating a divergence free initial condition and verifying div free

clear,clc %,close all

datapath = '/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-data/';
dirname = 'duct-knownmodes03';



NY = 65;
NZ = 49;

LY = 30;
LZ = 20;

filename = [datapath dirname '/grid.h5'];
gy = h5read(filename,'/grids/y');
gz = h5read(filename,'/grids/z');
gyf = (gy(1:end-1)+gy(2:end))/2;
gzf = (gz(1:end-1)+gz(2:end))/2;

% V = -cos(2*pi*gyf./LY *3/2)*sin(2*pi*gzf/LZ *3/2)'  *LY/10;
% W =  sin(2*pi*gyf./LY *3/2)*cos(2*pi*gzf/LZ *3/2)'  *LZ/10;
temp = 1./(cosh(gyf/LY*2*3)*cosh(gzf/LZ*2*3)');
V = -cos(2*pi*gyf./LY *3/2)*sin(2*pi*gzf/LZ *3/2)'  *LY/10 .*temp;
W =  sin(2*pi*gyf./LY *3/2)*cos(2*pi*gzf/LZ *3/2)'  *LZ/10 .*temp;



figure
imagesc(V')
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('V uniform')
% clim([-1 1])
% clim([-0.3 0.3])

figure
pcolor(gyf,gzf,V')
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('V nonuniform')


figure
imagesc(W')
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('W uniform')
% clim([-1 1])
% clim([-0.2 0.2])

figure
pcolor(gyf,gzf,W')
shading flat
colorbar
axis equal
xlabel('Y'), ylabel('Z')
title('W nonuniform')




% dudx = (U(2:end,:,:) - U(1:end-1,:,:)) ./ reshape((gx(2:end)-gx(1:end-1)),[length(gx)-1,1,1]);
dudx = 0;
dvdy = (V(2:end,:) - V(1:end-1,:)) ./ (gyf(2:end)-gyf(1:end-1));
dwdz = (W(:,2:end) - W(:,1:end-1)) ./ (gzf(2:end)-gzf(1:end-1))';

dvdyN = interp1(gzf,dvdy',gz(2:end-1))';
dwdzN = interp1(gyf,dwdz,gy(2:end-1));

div = dvdyN + dwdzN;
figure, imagesc(dvdyN), colorbar
figure, imagesc(dwdzN), colorbar
figure, imagesc(div), colorbar

% div = sum(dudx,'all') + sum(dvdy,'all') + sum(dwdz,'all');
% div = trapz(gzf,trapz(gy(2:end-1),abs(dvdy))) ...
%     + trapz(gyf,trapz(gz(2:end-1),abs(dwdz)'));
divsum = trapz(gz(2:end-1),trapz(gy(2:end-1),abs(div)));