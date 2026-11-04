% script for looking at various validation cases in duct
% right now for rectangular duct flow, pressure driven
% Poisseuille

clear,clc

datapath = '/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-data/';
dirname = 'duct-val-a1-16';

% NX = 16;
% NY = 65;
% NYP = 33;
% NZ = 49;
% NZP = 25;
% % NZ = 129;
% % NZP = 25;
% 
% NX = 32;
% NY = 129;
% % NYP = 33;
% NZ = 97;
% % NZP = 17;
% 
% NPROCZ = 2;
% NXP = NX/(2*NPROCZ);

filename = [datapath dirname '/grid.h5'];
gy = h5read(filename,'/grids/y');
gz = h5read(filename,'/grids/z');
gyf = (gy(1:end-1)+gy(2:end))/2;
gzf = (gz(1:end-1)+gz(2:end))/2;

file_info=h5info([datapath dirname '/mean.h5']);
att_info=file_info.Groups.Attributes;
nk=att_info.Value
% nk=31;

% LZ=20;
% NZ=96;
% gzf=linspace(0,LZ,NZ);


tii = [];
error = [];
mnmax = 500;


px = 1;
mu = 1;
a = gyf(end)/2;
b = gzf(end)/2;
u_a = zeros(length(gyf),length(gzf));
for j=1:length(gyf)
    for k=1:length(gzf)
        
        for m=0:mnmax
            for n=0:mnmax
                temp = (-1)^(m+n+1) /(2*m+1)/(2*n+1) ...
                    * (64/pi^4/mu)*(-1*px) ...
                    / ( ((2*m+1)/a)^2 + ((2*n+1)/b)^2 ) ...
                    * cos((2*m+1)*pi/2*(gyf(j)-a)/a) ...
                    * cos((2*n+1)*pi/2*(gzf(k)-b)/b);
                u_a(j,k) = u_a(j,k) + temp;
            end
        end

    end
end


filename=[datapath dirname '/movie.h5'];
figure
set(gcf,'Position',[0 300 1500 500])
set(gcf,'Color','w')
for k=1:nk
disp(k)


subplot(1,3,1)
varname=['/u_zy/' num2str(k,'%04d')];
A=h5read(filename,varname);
imagesc(gzf,gyf,A')
shading flat
axis tight
axis equal
colorbar
xlabel('Z'), ylabel('Y')
% u_zyprof = A(ceil(end/2),ceil(end/2));
u_xyprof = A(ceil(end/2),:);
u_xzprof = A(:,ceil(end/2));



subplot(1,3,2)
varname=['/u_xy/' num2str(k,'%04d')];
A=h5read(filename,varname);
imagesc(A')
shading flat
axis tight
axis equal
colorbar
xlabel('X'), ylabel('Y')
% u_xyprof = A(ceil(end/2),ceil(end/2));


subplot(1,3,3)
varname=['/u_xz/' num2str(k,'%04d')];
A=h5read(filename,varname);
imagesc(A')
shading flat
axis tight
axis equal
colorbar
xlabel('X'), ylabel('Z')
% u_xzprof = A(ceil(end/2),ceil(end/2));


% colormap(jet(256));
% colormap(gca,cmocean('balance'))

% clim([-2e-3 2e-3])
% clim([0 1])
% xlim([0 11])
% xlim([0 20])
% clim([-3 -1])
% axis

error(k) = max( max(abs(u_xyprof'-u_a(:,ceil(end/2)))),max(abs(u_xzprof'-u_a(ceil(end/2),:)))  );

t=h5read([datapath dirname '/mean.h5'],['/time/' num2str(k,'%04d')]);
text(0.1,0.8,['t=' num2str(t)],'units','normalized')
tii(k) = t;

M(k)=getframe(gcf);
clf;
end
close


%% final plots

% gyf = linspace(0,1,33);
% gzf = linspace(0,1,33);

% px = 1.0;
% mu = 1.0;
% a = gyf(end)/2;
% b = gzf(end)/2;
% 
% 
% u_a = zeros(length(gyf),length(gzf));
% for j=1:length(gyf)
%     for k=1:length(gzf)
%         
%         for m=0:mnmax
%             for n=0:mnmax
%                 temp = (-1)^(m+n+1) /(2*m+1)/(2*n+1) ...
%                     * (64/pi^4/mu)*(-1*px) ...
%                     / ( ((2*m+1)/a)^2 + ((2*n+1)/b)^2 ) ...
%                     * cos((2*m+1)*pi/2*(gyf(j)-a)/a) ...
%                     * cos((2*n+1)*pi/2*(gzf(k)-b)/b);
%                 u_a(j,k) = u_a(j,k) + temp;
%             end
%         end
% 
%     end
% end
figure
imagesc(gzf,gyf,u_a)
shading flat
colorbar


figure
subplot(1,2,1), hold on
plot(u_xyprof,gyf)
plot(u_a(:,ceil(end/2)),gyf)
xlabel('U'), ylabel('Y')
ylim([0 gyf(end)])

subplot(1,2,2), hold on
plot(u_xzprof,gzf)
plot(u_a(ceil(end/2),:),gzf)
xlabel('U'), ylabel('Z')
ylim([0 gzf(end)])


figure
subplot(1,2,1), hold on
plot(u_xyprof'-u_a(:,ceil(end/2)),gyf)
xlabel('U'), ylabel('Y')
ylim([0 gyf(end)])

subplot(1,2,2), hold on
plot(u_xzprof'-u_a(ceil(end/2),:),gzf)
xlabel('U'), ylabel('Z')
ylim([0 gzf(end)])


figure
semilogy(tii,error)

%%

gridvals = [16,32,64,128];
dzvals = 1./gridvals;
% errvals_a1 = [2.70e-4, 7.43e-5, 1.87e-5, 4.76e-6];
% errvals_a2 = [1.77e-4, 4.91e-5, 1.26e-5, 3.06e-6];
errvals_a1 = [2.70e-4, 7.41e-5, 1.86e-5, 4.64e-6];
errvals_a2 = [1.78e-4, 4.88e-5, 1.22e-5, 3.05e-6];

% x = 10:1:200;
% y = 0.1./(x.^2);
x = 1./(10:1:200);
y = 0.1 * (x.^2);
figure, hold on
plot(dzvals,errvals_a1,'LineWidth',1.5)
plot(dzvals,errvals_a2,'LineWidth',1.5)
scatter(dzvals,errvals_a1,50,[0 0 0],'filled','HandleVisibility', 'off')
scatter(dzvals,errvals_a2,50,[0 0 0],'filled','HandleVisibility', 'off')
plot(x,y,'--k','LineWidth',1.5)
set(gca,'YScale','log')
set(gca,'XScale','log')
legend("$a/b=1$","$a/b=2$","$\sim \Delta z^2$")
xlabel('$\Delta z$')
ylabel('error')
legend('Location','NorthWest')



% FS = 10;      % Font Size
% FS_l = 7;     % Font Size - legends
% FN = 'Times New Roman';  % Font Name
% picturewidth1 = 6.49; % inches
% picturewidth2 = 3.5; % inches
% picturewidth3 = 5.0; % inches
% hw_ratio = 1.0;
set(findall(gcf,'-property','FontSize'),'FontSize',10)
% set(findall(gcf,'-property','Box'),'Box','on')
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex') 
set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(gcf,'Units','inches','Position',[1 1 3.5 4])
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','inches','PaperSize',[pos(3), pos(4)])
savepath='/Users/vincent/Documents/_School/Berkeley/Research/Thesis/CH3/Figures/';
print(gcf,[savepath 'poiseuille'],'-dpdf','-vector','-fillpage');

% return
% %% Code to save movie (optional)
% filename_mov=[datapath dirname '/movie_01_wide_zy.mp4'];
% clear v
% v = VideoWriter(filename_mov, 'MPEG-4');
% v.FrameRate = 10;
% open(v)
% writeVideo(v,M)
% close(v)
