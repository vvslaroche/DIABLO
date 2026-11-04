% copy of debug7 but for comparing 2 movies together

clear,clc

datapath = '/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-data/';
dirname = 'duct-tests';


NX = 16;
NY = 65;
NYP = 33;
NZ = 49;
NZP = 25;

NPROCZ = 2;
NXP = NX/(2*NPROCZ);

filename = [datapath 'duct-knownmodes03-ai2' '/grid.h5'];
gy = h5read(filename,'/grids/y');
gz = h5read(filename,'/grids/z');
gyf = (gy(1:end-1)+gy(2:end))/2;
gzf = (gz(1:end-1)+gz(2:end))/2;

nk=201;

% LZ=4;
% 
% NZ=16;
% z=linspace(0,LZ,NZ);

filename1=[datapath dirname '/05-movie.h5'];
filename2=[datapath dirname '/04-movie.h5'];
figure
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

varname=['/v_zy/' timename];
A=h5read(filename1,varname);
B=h5read(filename2,varname);
imagesc(A-B)
% pcolor(z,gyf,A')
% caxis([-1.5 1.5]);
% axis tight
% shading interp
shading flat
xlabel('Y'), ylabel('Z')
% colormap(jet(256));
colorbar
% clim([-1e-3 1e-3])
% clim([-1.5 1.5])
M(k)=getframe(gcf);
clf;
end
close

return
%% Code to save movie (optional)
filename_mov=[datapath dirname '/movie_zy.mp4'];
clear v
v = VideoWriter(filename_mov, 'MPEG-4');
v.FrameRate = 10;
open(v)
writeVideo(v,M)
close(v)
