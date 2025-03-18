% This script shows how to load in 2D slices and make a movie of the simulation output
% Run after readmean.m
LX=13.98;
NX=128;

x=linspace(0,LX,NX);

filename=[base_dir '/movie.h5'];

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

varname=['/th04_xy/' timename];
% varname=['/u_xy/' timename];
%varname=['/nu_t_xy/' timename];

A=h5read(filename,varname);

pcolor(x,gyf,A'); shading interp;
%caxis([-1.5 1.5]);


axis tight
ylim([-10 10])
shading interp
colormap(jet(256));
colorbar
M(k)=getframe(gcf);
clf;
end
close


% saving movie code
filename_mov=[base_dir '/movie.mp4'];
clear v
v = VideoWriter(filename_mov, 'MPEG-4');
v.FrameRate = 20;
open(v)
writeVideo(v,M)
close(v)

