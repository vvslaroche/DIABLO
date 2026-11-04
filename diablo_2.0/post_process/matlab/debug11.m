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