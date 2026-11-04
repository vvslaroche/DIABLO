% code to calculate cumulative mixing efficiency eta_c
% from 2D movie slices
% intended for use on duct geometry gravity current simulations
% for LANL visit

clear, clc


datapath = '/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-data/';
dirname = 'duct-GC-NZ4096-NU5.0e-5';


% 
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


filename = [datapath dirname '/grid.h5'];
gy = h5read(filename,'/grids/y');
gz = h5read(filename,'/grids/z');
gyf = (gy(1:end-1)+gy(2:end))/2;
gzf = (gz(1:end-1)+gz(2:end))/2;

dy = gyf(5)-gyf(4);
dz = gzf(5)-gzf(4);

nk=133;

Lz = 1;
Ri_b = 1;
Pr = 1;
Re = 1/(5.0e-5);
normalization_factor = length(gzf)*length(gyf);


tii = [];

nbins = 101;
buoypdf = zeros(nbins,nk);
bincenter = linspace(0,1,nbins);
binwidth = 1/(nbins-1);
bincount = zeros(nbins,nk);

filename=[datapath dirname '/movie.h5'];
figure
% set(gcf,'Position',[0 100 1500 180])
% set(gcf,'Color','w')
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
A=h5read(filename,varname);

for J=1:length(gyf)
    for K=1:length(gzf)
        for b=1:nbins
            if A(K,J)>bincenter(b)-binwidth/2 && A(K,J)<bincenter(b)+binwidth/2
                bincount(b,k) = bincount(b,k) +1;
            end
        end
        if A(K,J)<bincenter(1)-binwidth/2
            bincount(1,k) = bincount(1,k) +1;
        end
        if A(K,J)>bincenter(end)+binwidth/2
            bincount(end,k) = bincount(end,k) +1;
        end
    end
end

% check that all values are being binned
if sum(bincount(:,k))~=length(gyf)*length(gzf)
    disp('uncounted values detected')
    return
end

t=h5read([datapath dirname '/mean.h5'],['/time/' timename]);
tii(k) = t;

thme=h5read([datapath dirname '/mean.h5'],['/thme01/' timename]);



Z=Lz*cumtrapz(bincenter,bincount(:,k)/normalization_factor); 

% calculate total PE, BPE, and APE
PE_b(k)=-Ri_b*trapz(Z,(bincenter'.*Z))/Lz;
D_pt(k)=-Ri_b*(-1-1)/(Re*Pr*Lz)*tii(k); % say for now delta(t) is approximately 5e-2
% D_p(ts)=-Ri_b*(-1-1)/(Re*Pr*Lz);
PE(k)=-Ri_b*trapz(gyf,gyf.*thme)/Lz;


% scatter(bincenter,bincount(:,k))



% M(k)=getframe(gcf);
% clf;
end
close


% eta = (M+e)/M
% M = phi_d - phi_i
% phi_d = d/dt of background PE
% phi_i = increase of background PE due to diffusion


%% energetics
PE_a = PE - PE_b;

epsilon_tot = epsilon + 1/Re*(dudy.^2 + dwdy.^2);
epsilon_tot_int=squeeze(trapz(gy,epsilon_tot,1));

m = diff(movmean(PE_b - D_pt,7));
figure
plot(tii(2:end),m)

figure
plot(tii,epsilon_tot_int/30,tiif,m)
legend("epsilon","m")

eta = m ./ (m + epsilon_tot_int(2:end)/30);
figure; plot(tiif,movmean(eta,7))
title('mixing efficiency')

cumeta = trapz(tiif(20:end),m(20:end)) / trapz(tiif(20:end),m(20:end) + epsilon_tot_int(21:end));

figure, hold on
plot(tii,PE)
plot(tii,PE_a)
plot(tii,PE_b)


%% plotting
figure, hold on
plot(bincenter,cumtrapz(bincenter,bincount(:,1)/normalization_factor))
plot(bincenter,cumtrapz(bincenter,bincount(:,50)/normalization_factor))
plot(bincenter,cumtrapz(bincenter,bincount(:,100)/normalization_factor))
plot(bincenter,cumtrapz(bincenter,bincount(:,133)/normalization_factor))



return
%% Code to save movie (optional)
filename_mov=[datapath dirname '/movie_zy.mp4'];
clear v
v = VideoWriter(filename_mov, 'MPEG-4');
v.FrameRate = 10;
open(v)
writeVideo(v,M)
close(v)

