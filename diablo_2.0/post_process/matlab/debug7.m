% basically plot_movie_zy but for debugging purposes

clear,clc

datapath = '/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-data/';
% dirname = 'duct-knownmodes03-ai2';
% dirname = 'duct-GC6-NZ2048-hpc-NUe-5';
% dirname = 'duct-GC-NZ8192-NU1.0e-5-savio_02';
% dirname = 'duct-GC-NZ6144-NU2.5e-5';
dirname = 'duct-GC-NZ4096-NU5.0e-5';
% dirname = 'duct-GC-NZ2048-NU1.0e-4-DA0.10';
% dirname = 'duct-KH-NY128-NU1.0e-3-xme';


NX = 16;
NY = 65;
NYP = 33;
NZ = 49;
NZP = 25;
% NZ = 129;
% NZP = 25;

NX = 32;
NY = 129;
% NYP = 33;
NZ = 97;
% NZP = 17;

NPROCZ = 2;
NXP = NX/(2*NPROCZ);

filename = [datapath dirname '/grid.h5'];
gy = h5read(filename,'/grids/y');
gz = h5read(filename,'/grids/z');
gyf = (gy(1:end-1)+gy(2:end))/2;
gzf = (gz(1:end-1)+gz(2:end))/2;

nk=501;

% LZ=20;
% NZ=96;
% gzf=linspace(0,LZ,NZ);

bottom_row = [];
tii = [];

maxE = 0;

filename=[datapath dirname '/movie.h5'];
figure
set(gcf,'Position',[0 100 1500 180])
set(gcf,'Color','w')
for k=1:nk
k
  if (k<10)
    timename=['000' int2str(k)];
  elseif (k<100) 
    timename=['00' int2str(k)];
  elseif (k<1000)
    timename=['0' int2str(k)];
  else
    timename=[int2str(k)];
  end

% varname=['/epsilon_zy/' timename];
varname=['/th01_zy/' timename];
A=h5read(filename,varname);

% B=h5read(filename,varname);
% imagesc(A)
% A(A<0)=1;
% pcolor(gzf,gyf,A'-(1-B'))
% pcolor(gzf,gyf,log10(A)')
pcolor(gzf,gyf,A')
% pcolor(gzf,gyf,fliplr(A'))
shading flat

axis tight
axis equal

% if max(max(A))>maxE
%     maxE = max(max(A))
% end

% xlabel('Y'), ylabel('Z')
xlabel('Z'), ylabel('Y')
% colormap(jet(256));
colormap(gca,cmocean('balance'))
colorbar
% clim([-2e-3 2e-3])
clim([0 1])
% xlim([0 11])
xlim([0 20])
% clim([-3 -1])

t=h5read([datapath dirname '/mean.h5'],['/time/' timename]);
text(0.1,0.8,['t=' num2str(t)],'units','normalized','color','white')
tii(k) = t;
bottom_row = [bottom_row; A(:,2)'];

M(k)=getframe(gcf);
clf;
end
close

return
%% Code to save movie (optional)
filename_mov=[datapath dirname '/movie_01_wide_zy.mp4'];
clear v
v = VideoWriter(filename_mov, 'MPEG-4');
v.FrameRate = 10;
open(v)
writeVideo(v,M)
close(v)


%% code to calculate Kolmogorov scale and required grid sizing (optional)
maxE = 0.01; % approximate max epsilon from movie
NU = 5.0e-5;
% figure; set(gcf,'defaultaxescolororder',parula(8)); plot(epsilon(:,51:50:end),gyf)
% maxE = max(max(epsilon(:,20:end)));
L_k = (NU^3/maxE)^(1/4)
% L_b = L_k/sqrt(Pr)      % batchelor scale
% needed_points_Y = ceil(30/(2.5*L_b))
% needed_points_X_H = ceil(15.22/(2.5*L_b))
% needed_points_X_KH = ceil(13.98/(2.5*L_b))
% needed_points_Z = ceil(5/(2.5*L_b))
% needed_points_Y = ceil(1/(2.5*L_k))
% needed_points_X_H = ceil(0.5/(2.5*L_b))
% needed_points_X_KH = ceil(13.98/(2.5*L_b))
needed_points_Z = ceil(20/(2.5*L_k))
% needed_points_Z = ceil(20/(5*L_k))


%% code to calculate gravity current speed (optional)
figure
pcolor(tii,gzf,bottom_row')
shading flat
colorbar
xlabel('time')
ylabel('z location of front')
clim([0 1])


% test = find(bottom_row<0.99,1,'first');
for i = 1:size(bottom_row,1)
    loc(i) = gzf(find(bottom_row(i,:) < 0.99, 1, 'last'));
end
figure, hold on
plot(tii,loc)

loci1 = find(tii>10,1,'first'); %25
loci2 = find(tii>20,1,'first'); %50
h=fit(tii(loci1:loci2)',loc(loci1:loci2)','poly1');
plot(tii(loci1:loci2),h.p1*tii(loci1:loci2)+h.p2)
h.p1