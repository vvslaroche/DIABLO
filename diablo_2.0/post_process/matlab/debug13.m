% basically plot_movie_zy but for calculation of 
% x averaged TKE budget

clear,clc

datapath = '/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-data/';
dirname = 'duct-KH-NY128-NU1.0e-3-xme';


filename = [datapath dirname '/grid.h5'];
gy = h5read(filename,'/grids/y');
gz = h5read(filename,'/grids/z');
gyf = (gy(1:end-1)+gy(2:end))/2;
gzf = (gz(1:end-1)+gz(2:end))/2;
LY = gyf(end)-gyf(1);
LZ = gzf(end)-gzf(1);

RI = 0.15;

nk=351;



tii = zeros(1,nk);
tke = zeros(1,nk);
shear_prod = zeros(1,nk);
buoy_prod = zeros(1,nk);
epsilon = zeros(1,nk);

filename=[datapath dirname '/movie.h5'];
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

% varname=['/thxme01_zy/' timename];
% A=h5read(filename,varname);
% % imagesc(A)
% % A(A<0)=1;
% pcolor(gzf,gyf,A')
% % pcolor(gzf,gyf,fliplr(A'))
% shading flat

temp = 0;
varname=['/urmsxme_zy/' timename];
A=h5read(filename,varname);
temp = temp + 0.5*(A.^2);
varname=['/vrmsxme_zy/' timename];
A=h5read(filename,varname);
temp = temp + 0.5*(A.^2);
varname=['/wrmsxme_zy/' timename];
A=h5read(filename,varname);
temp = temp + 0.5*(A.^2);
tke(k) = trapz(gyf,trapz(gzf,temp))/LY/LZ;


temp = 0;
varname=['/uvxme_zy/' timename];
A=h5read(filename,varname);
varname=['/duxmedy_zy/' timename];
B=h5read(filename,varname);
temp = temp + A.*B;
varname=['/uwxme_zy/' timename];
A=h5read(filename,varname);
varname=['/duxmedz_zy/' timename];
B=h5read(filename,varname);
temp = temp + A.*B;
varname=['/vrmsxme_zy/' timename];
A=h5read(filename,varname);
varname=['/dvxmedy_zy/' timename];
B=h5read(filename,varname);
temp = temp + (A.^2).*B;
varname=['/wvxme_zy/' timename];
A=h5read(filename,varname);
varname=['/dvxmedz_zy/' timename];
B=h5read(filename,varname);
temp = temp + A.*B;
varname=['/wvxme_zy/' timename];
A=h5read(filename,varname);
varname=['/dwxmedy_zy/' timename];
B=h5read(filename,varname);
temp = temp + A.*B;
varname=['/wrmsxme_zy/' timename];
A=h5read(filename,varname);
varname=['/dwxmedz_zy/' timename];
B=h5read(filename,varname);
temp = temp + (A.^2).*B;
temp = -temp;
shear_prod(k) = trapz(gyf,trapz(gzf,temp))/LY/LZ;


varname=['/thvxme01_zy/' timename];
A=h5read(filename,varname);
temp = A;
buoy_prod(k) = -RI(1) * trapz(gyf,trapz(gzf,temp))/LY/LZ;


varname=['/epsxme_zy/' timename];
A=h5read(filename,varname);
temp = A;
epsilon(k) = trapz(gyf,trapz(gzf,temp))/LY/LZ;

% prod_xme =
%  -( uv*dUdy + uw*dUdz
%   + vv*dVdy + vw*dVdz
%   + vw*dWdy + ww*dWdz )
% axis tight
% axis equal
% 
% % xlabel('Y'), ylabel('Z')
% xlabel('Z'), ylabel('Y')
% % colormap(jet(256));
% % colormap(gca,cmocean('algae'))
% colorbar
% % clim([-2e-3 2e-3])
% clim([0 1])
% % xlim([0 11])
% 
tii(k)=h5read([datapath dirname '/mean.h5'],['/time/' timename]);
% text(0.1,0.8,['t=' num2str(t)],'units','normalized')
% tii(k) = t;
% bottom_row = [bottom_row; A(:,2)'];
% 
% M(k)=getframe(gcf);
% clf;
end
% close

dtkedt = diff(tke)./diff(tii);
tiif = (tii(2:end)+tii(1:end-1))/2;


figure
hold on
plot(tiif,dtkedt, 'LineWidth',1)
plot(tii,shear_prod, 'LineWidth',1)
plot(tii,-epsilon, 'LineWidth',1)
plot(tii,-buoy_prod, 'LineWidth',1)
plot(tii,shear_prod - epsilon - buoy_prod,'--k','LineWidth',2)
legend("d(tke)/dt","shear prod","- epsilon","- buoy prod","summed",'FontSize',20)
xlabel('time')

return
%% Code to save movie (optional)
filename_mov=[datapath dirname '/movie_zy.mp4'];
clear v
v = VideoWriter(filename_mov, 'MPEG-4');
v.FrameRate = 10;
open(v)
writeVideo(v,M)
close(v)


%%
% check TKE budget for DUCT
% run after readmean

%     tke(j,k)=0.5*(urms(j,k)^2.+vrms(j,k)^2.+wrms(j,k)^2.);
for k=1:nk
  for j=1:NY
    tke(j,k)=0.5*(urms(j,k)^2.+vrms(j,k)^2.+wrms(j,k)^2.);
    shear_prod(j,k)= -(uv(j,k)*dudy(j,k) + wv(j,k)*dwdy(j,k)); 
%         - (urms.*dudx + vrms.*dvdy + wrms.*dwdz + uv.*dvdx + uw.*dudz + uw.*dwdx + wv.*dvdz)
    buoy_prod(j,k)=-RI(1) * thv(j,k,1);
  end
end
% shear_prod_int=trapz(gyf,shear_prod,1);
shear_prod_int=trapz(gyf(1:end-1),shear_prod(1:end-1,:),1) /(gyf(end-1)-gyf(1));
epsilon_int=trapz(gyf,epsilon,1) /(gyf(end)-gyf(1));
buoy_prod_int=trapz(gyf,buoy_prod,1) /(gyf(end)-gyf(1));

tke_int=trapz(gyf,tke,1);

dtkeintdt = diff(tke_int/(gyf(end)-gyf(1)))./diff(tii);
tiif = (tii(2:end)+tii(1:end-1))/2;


figure
hold on
plot(tiif,dtkeintdt, 'LineWidth',1)
% plot(tii,shear_prod_int, tii,-epsilon_int, tii,buoy_prod_int)
% plot(tii,shear_prod_int - epsilon_int - buoy_prod_int,'--','LineWidth',2)
plot(tii,shear_prod_int, tii,-epsilon_int, tii,-buoy_prod_int, 'LineWidth',1)
plot(tii,shear_prod_int - epsilon_int - buoy_prod_int,'--k','LineWidth',2)
legend("d(tke)/dt","shear prod","- epsilon","- buoy prod","summed",'FontSize',20)
% xlim([tii(10) tii(end)/3600])
xlabel('time')
% exportgraphics(gcf,'tke_budget.png')