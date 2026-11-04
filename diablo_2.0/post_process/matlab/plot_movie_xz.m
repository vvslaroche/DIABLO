% This script shows how to load in 2D slices and make a movie
% Run after readmean_h5.m
LX=0.5;
NX=144;
x=linspace(0,LX,NX);
LZ=20;
NZ=6145;
z=linspace(0,LZ,NZ);

filename=[base_dir '/movie.h5'];


figure
set(gcf,'Position',[30 100 1450 100])
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

varname=['/th02_xz/' timename];
% varname=['/u_xz/' timename];
A=h5read(filename,varname);

% pcolor(x,z,A')
pcolor(z,x,A)
% caxis([-1.5 1.5]);
clim([0 1])
axis tight
shading interp
xlabel('X'), ylabel('Z')
% colormap(jet(256));
colormap(cmocean('algae'))
colorbar
% t=h5read([datapath dirname '/mean.h5'],['/time/' timename]);
text(0.1,0.8,['t=' num2str(tii(k))],'units','normalized','Color','w')


M(k)=getframe(gcf);
clf;
end
close


return
%% Code to save movie (optional)
filename_mov=[base_dir '/movie_xz.mp4'];
clear v
v = VideoWriter(filename_mov, 'MPEG-4');
v.FrameRate = 10;
open(v)
writeVideo(v,M)
close(v)
