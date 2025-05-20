% This script shows how to load in 2D slices and make a movie
% Run after readmean_h5.m
LZ=4;
NZ=17;
z=linspace(0,LZ,NZ);

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

varname=['/th01_zy/' timename];
varname=['/u_zy/' timename];
A=h5read(filename,varname);

pcolor(z,gyf,A')
% caxis([-1.5 1.5]);
axis tight
shading interp
xlabel('Z'), ylabel('Y')
colormap(jet(256));
colorbar
M(k)=getframe(gcf);
clf;
end
close


%% Code to save movie (optional)
filename_mov=[base_dir '/movie_zy.mp4'];
clear v
v = VideoWriter(filename_mov, 'MPEG-4');
v.FrameRate = 20;
open(v)
writeVideo(v,M)
close(v)
