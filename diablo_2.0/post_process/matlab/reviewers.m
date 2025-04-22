
cmapfile = importdata('/Users/vincent/Documents/colorbars/MPL_PuOr.rgb',' ',2);
cmap1 = cmapfile.data;
cmap1 = [cmap1(1:64,:); [1 1 1]; cmap1(65:128,:)];
cmap1 = flipud(cmap1);
cmapfile = importdata('/Users/vincent/Documents/colorbars/cmocean_algae.rgb',' ',2);
cmap2 = cmapfile.data/255;
cmapfile = importdata('/Users/vincent/Documents/colorbars/NCV_blu_red.rgb',' ',2);
cmap3 = cmapfile.data/255;

c2 = [0.5 0.75 0.9 0.99];
c2 = [1-flip(c2(2:end)) c2];
c1 = 2*c2-1;

gyf2 = gyf-15;

figure
subplot(2,1,1)
pcolor(tii,gyf2,thv(:,:,1)); shading flat;
cb=colorbar; colormap(gca,cmap1); clim([-5e-3 5e-3]);
hold on; contour(tii,gyf2,thme(:,:,1),c1,'k-','LineWidth',.2);
contour(tii,gyf2,thme(:,:,1),[0.0 0.0],'k-','LineWidth',1.5);
ylabel(cb,'$\phi_a$','Interpreter','Latex')
ylabel('z'); xlabel('t');
axis([0 300 -19 -11])

subplot(2,1,2)
pcolor(tii,gyf2,thv(:,:,2)); shading flat;
cb=colorbar; colormap(gca,cmap1); clim([-5e-3 5e-3]);
hold on; contour(tii,gyf2,thme(:,:,2),c2,'k-','LineWidth',.2);
contour(tii,gyf2,thme(:,:,2),[0.5 0.5],'k-','LineWidth',1.5);
ylabel(cb,'$\phi_a$','Interpreter','Latex')
ylabel('z'); xlabel('t');
axis([0 300 -19 -11])



c2 = [0.0100    0.1000    0.2500      0.7500    0.9000    0.9900];

figure
subplot(4,1,1)
pcolor(tii,gyf2,thme(:,:,1)); shading flat
colorbar, colormap(gca,cmap3), clim([-1 1])
hold on; contour(tii,gyf2,thme(:,:,1),c1,'k-','LineWidth',.2);
contour(tii,gyf2,thme(:,:,1),[0.0 0.0],'k-','LineWidth',1.5);
ylabel('z'); xlabel('t'); title('Buoyancy')
axis([0 400 -19 -11])

subplot(4,1,2)
pcolor(tii,gyf2,thme(:,:,2)); shading flat
colorbar, colormap(gca,cmap2), clim([0 1])
hold on; contour(tii,gyf2,thme(:,:,2),c2,'k-','LineWidth',.2);
contour(tii,gyf2,thme(:,:,4),[0.5 0.5],'k-','LineWidth',1.5);
ylabel('z'); xlabel('t'); title('Phytoplankton, Da=0.1')
axis([0 400 -19 -11])

subplot(4,1,3)
pcolor(tii,gyf2,thme(:,:,4)); shading flat
colorbar, colormap(gca,cmap2), clim([0 1])
hold on; contour(tii,gyf2,thme(:,:,4),c2,'k-','LineWidth',.2);
contour(tii,gyf2,thme(:,:,4),[0.5 0.5],'k-','LineWidth',1.5);
ylabel('z'); xlabel('t'); title('Phytoplankton, Da=1')
axis([0 400 -19 -11])

subplot(4,1,4)
pcolor(tii,gyf2,thme(:,:,6)); shading flat
colorbar, colormap(gca,cmap2), clim([0 1])
hold on; contour(tii,gyf2,thme(:,:,6),c2,'k-','LineWidth',.2);
contour(tii,gyf2,thme(:,:,4),[0.5 0.5],'k-','LineWidth',1.5);
ylabel('z'); xlabel('t'); title('Phytoplankton, Da=10')
axis([0 400 -19 -11])


set(gcf,'Position',[100 100 400 800])



%% extra plots of P'w'
figure
subplot(4,1,1)
pcolor(tii,gyf2,thv(:,:,1)); shading flat
colorbar, colormap(gca,cmap1), clim([-5e-3 5e-3]);
hold on; contour(tii,gyf2,thme(:,:,1),c1,'k-','LineWidth',.2);
contour(tii,gyf2,thme(:,:,1),[0.0 0.0],'k-','LineWidth',1.5);
ylabel('z'); xlabel('t'); title('Buoyancy')
axis([0 400 -19 -11])

subplot(4,1,2)
pcolor(tii,gyf2,thv(:,:,2)); shading flat
colorbar, colormap(gca,cmap1), clim([-5e-3 5e-3]);
hold on; contour(tii,gyf2,thme(:,:,2),c2,'k-','LineWidth',.2);
contour(tii,gyf2,thme(:,:,4),[0.5 0.5],'k-','LineWidth',1.5);
ylabel('z'); xlabel('t'); title('Phytoplankton, Da=0.1')
axis([0 400 -19 -11])

subplot(4,1,3)
pcolor(tii,gyf2,thv(:,:,4)); shading flat
colorbar, colormap(gca,cmap1), clim([-5e-3 5e-3]);
hold on; contour(tii,gyf2,thme(:,:,4),c2,'k-','LineWidth',.2);
contour(tii,gyf2,thme(:,:,4),[0.5 0.5],'k-','LineWidth',1.5);
ylabel('z'); xlabel('t'); title('Phytoplankton, Da=1')
axis([0 400 -19 -11])

subplot(4,1,4)
pcolor(tii,gyf2,thv(:,:,6)); shading flat
colorbar, colormap(gca,cmap1), clim([-5e-3 5e-3]);
hold on; contour(tii,gyf2,thme(:,:,6),c2,'k-','LineWidth',.2);
contour(tii,gyf2,thme(:,:,4),[0.5 0.5],'k-','LineWidth',1.5);
ylabel('z'); xlabel('t'); title('Phytoplankton, Da=10')
axis([0 400 -19 -11])


set(gcf,'Position',[100 100 400 800])