% This script shows how to load in 2D slices and make a movie
% Run after readmean_h5.m
LX=14;
NX=384;
x=linspace(0,LX,NX);

filename=[base_dir '/movie.h5'];


figure
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

varname=['/th01_xy/' timename];
% varname=['/u_xy/' timename];
A=h5read(filename,varname);
% B=h5read('../../../../DIABLO-data/duct-chan-knownmodes-copy2/movie.h5',varname);

% pcolor(x,gyf,A'-B')

pcolor(x,gyf,A')
% imagesc(A')
% caxis([-1.5 1.5]);
axis tight
shading interp
xlabel('X'), ylabel('Y')
colormap(jet(256));
colorbar
ylim([-5 5])

text(0.1,0.8,['t=' num2str(tii(k))],'units','normalized')

M(k)=getframe(gcf);
clf;
end
close


return
%% Code to save movie (optional)
filename_mov=[base_dir '/movie_xy.mp4'];
clear v
v = VideoWriter(filename_mov, 'MPEG-4');
v.FrameRate = 20;
open(v)
writeVideo(v,M)
close(v)
