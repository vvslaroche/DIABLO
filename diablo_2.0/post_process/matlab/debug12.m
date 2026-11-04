% calculate TKE budget from 3D saves
% should be flexible enough to run for either CHAN or DUCT

clear,clc

datapath = '/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-data/';
dirname = 'duct-KH-hpc-highres';

NUM_PER_DIR = 2;
SAVE_FLOW_INT = 200;
NX = 96;
NY = 129;
NZ = 65;
LX = 14;
LZ = 10;
RI = 0.15;
NU = 0.001;

% end user input

gx = linspace(0,LX,NX)';
gxf = (gx(2:end)+gx(1:end-1))/2;

gz = linspace(0,LZ,NZ)';
gzf = (gz(2:end)+gz(1:end-1))/2;

filename = [datapath dirname '/grid.h5'];
gy = h5read(filename,'/grids/y');
gyf = (gy(1:end-1)+gy(2:end))/2;
LY = gyf(end)-gyf(1);


if NUM_PER_DIR==1
    gz = h5read(filename,'/grids/z');
    gzf = (gz(1:end-1)+gz(2:end))/2;
%     LZ = gzf(end)-gzf(1);
    gz = gzf;
end


time = nan(200,1);
tke = nan(200,1);
shear_prod = nan(200,1);
buoy_prod = nan(200,1);
epsilon = nan(200,1);
urms_3D = nan(200,129);
for i=1:200
    filename=[datapath dirname '/output/out' num2str(SAVE_FLOW_INT*i,'%05.0f') '.h5'];
    disp(['Reading out' num2str(SAVE_FLOW_INT*i,'%05.0f') '.h5'])
    try
        time(i) = h5readatt(filename,'/Timestep','Time');
    catch
        disp(['out' num2str(SAVE_FLOW_INT*i,'%05.0f') '.h5 does not exist. Exiting...'])
        break
    end
    U = h5read(filename,'/Timestep/U');
    V = h5read(filename,'/Timestep/V');
    W = h5read(filename,'/Timestep/W');
    
    

    % in CHAN, all variables on GYF, (and GX,GZ)
    % in DUCT, all variables on GYF,GZF (and GX)

    % calculate means
    ume = trapz(gz,trapz(gx,U,1),3) /LX/LZ;
%     ume = trapz(gzf,trapz(gx,U,1),3);
    vme = trapz(gz,trapz(gx,V,1),3) /LX/LZ;
    wme = trapz(gz,trapz(gx,W,1),3) /LX/LZ;
    

    % calculate turbulent fluctuations
%     uu = (U-ume).^2;
%     vv = (V-vme).^2;
%     ww = (W-wme).^2;
%     uv = (U-ume).*(V-vme);
%     wv = (W-wme).*(V-vme);
    uu = trapz(gz,trapz(gx, (U-ume).^2 ,1),3) /LX/LZ;
    vv = trapz(gz,trapz(gx, (V-vme).^2 ,1),3) /LX/LZ;
    ww = trapz(gz,trapz(gx, (W-wme).^2 ,1),3) /LX/LZ;
    uv = trapz(gz,trapz(gx, (U-ume).*(V-vme) ,1),3) /LX/LZ;
    wv = trapz(gz,trapz(gx, (W-wme).*(V-vme) ,1),3) /LX/LZ;

    urms_3D(i,:) = sqrt(uu);

    % calculate gradients
    dudy = (ume(2:end)-ume(1:end-1))./(gyf(2:end)-gyf(1:end-1))';
    dwdy = (wme(2:end)-wme(1:end-1))./(gyf(2:end)-gyf(1:end-1))';

%     tke(i) = trapz(gz, ...
%              trapz(gyf, ...
%              trapz(gx,uu + vv + ww))) /LX/LY/LZ;
    tke(i) = trapz(gyf,uu + vv + ww) /LY;

    shear_prod(i) = - trapz((gyf(2:end)+gyf(1:end-1))/2, ...
              (uv(2:end)+uv(1:end-1))/2 .* dudy ...
             +(wv(2:end)+wv(1:end-1))/2 .* dwdy ) /LY;


%     clear U V W uu vv ww
    % buoy prod term
    TH = h5read(filename,'/Timestep/TH01');
    thme = trapz(gz,trapz(gx,TH,1),3) /LX/LZ;
    thv = trapz(gz,trapz(gx, (TH-thme).*(V-vme) ,1),3) /LX/LZ;
    buoy_prod(i) = trapz(gyf,-RI*thv) /LY;


    % epsilon
    epsilon(i) = 0;
    eps3D = 0;
    % d/dx terms
    temp = ((U(2:end,:,:)-U(1:end-1,:,:))./(gx(2:end)-gx(1:end-1))).^2;
    temp = (temp(:,2:end,:)+temp(:,1:end-1,:))/2;
    temp = (temp(:,:,2:end)+temp(:,:,1:end-1))/2;
    eps3D = eps3D + temp;
    temp = ((V(2:end,:,:)-V(1:end-1,:,:))./(gx(2:end)-gx(1:end-1))).^2;
    temp = (temp(:,2:end,:)+temp(:,1:end-1,:))/2;
    temp = (temp(:,:,2:end)+temp(:,:,1:end-1))/2;
    eps3D = eps3D + temp;
    temp = ((W(2:end,:,:)-W(1:end-1,:,:))./(gx(2:end)-gx(1:end-1))).^2;
    temp = (temp(:,2:end,:)+temp(:,1:end-1,:))/2;
    temp = (temp(:,:,2:end)+temp(:,:,1:end-1))/2;
    eps3D = eps3D + temp;
    % d/dz terms
    gzn = permute(gz,[3,2,1]);
    temp = ((U(:,:,2:end)-U(:,:,1:end-1))./(gzn(2:end)-gzn(1:end-1))).^2;
    temp = (temp(:,2:end,:)+temp(:,1:end-1,:))/2;
    temp = (temp(2:end,:,:)+temp(1:end-1,:,:))/2;
    eps3D = eps3D + temp;
    temp = ((V(:,:,2:end)-V(:,:,1:end-1))./(gzn(2:end)-gzn(1:end-1))).^2;
    temp = (temp(:,2:end,:)+temp(:,1:end-1,:))/2;
    temp = (temp(2:end,:,:)+temp(1:end-1,:,:))/2;
    eps3D = eps3D + temp;
    temp = ((W(:,:,2:end)-W(:,:,1:end-1))./(gzn(2:end)-gzn(1:end-1))).^2;
    temp = (temp(:,2:end,:)+temp(:,1:end-1,:))/2;
    temp = (temp(2:end,:,:)+temp(1:end-1,:,:))/2;
    eps3D = eps3D + temp;
    % d/dy terms
    u = U-ume;
    temp = ((u(:,2:end,:)-u(:,1:end-1,:))./(gyf(2:end)-gyf(1:end-1))').^2;
    temp = (temp(2:end,:,:)+temp(1:end-1,:,:))/2;
    temp = (temp(:,:,2:end)+temp(:,:,1:end-1))/2;
    eps3D = eps3D + temp;
    v = V-vme;
    temp = ((v(:,2:end,:)-v(:,1:end-1,:))./(gyf(2:end)-gyf(1:end-1))').^2;
    temp = (temp(2:end,:,:)+temp(1:end-1,:,:))/2;
    temp = (temp(:,:,2:end)+temp(:,:,1:end-1))/2;
    eps3D = eps3D + temp;
    w = W-wme;
    temp = ((w(:,2:end,:)-w(:,1:end-1,:))./(gyf(2:end)-gyf(1:end-1))').^2;
    temp = (temp(2:end,:,:)+temp(1:end-1,:,:))/2;
    temp = (temp(:,:,2:end)+temp(:,:,1:end-1))/2;
    eps3D = eps3D + temp;

    gz_temp = (gz(2:end)+ gz(1:end-1))/2;
    gyf_temp = (gyf(2:end)+gyf(1:end-1))/2;
    gx_temp = (gx(2:end)+ gx(1:end-1))/2;
    epsilon(i) = NU* trapz(gz_temp, ...
                     trapz(gyf_temp, ...
                     trapz(gx_temp,eps3D))) ...
       /(gx_temp(end)-gx_temp(1)) /(gyf_temp(end)-gyf_temp(1)) /(gz_temp(end)-gz_temp(1));

end


% return
%% plotting
figure, hold on
plot((time(2:end)+time(1:end-1))/2,(tke(2:end)-tke(1:end-1))./(time(2:end)-time(1:end-1)))
% plot(tii,shear_prod, tii,-epsilon_int, tii,-buoy_prod_int)
% plot(tii,shear_prod - epsilon_int - buoy_prod_int,'--k','LineWidth',2)
plot(time,shear_prod)
plot(time,-epsilon)
plot(time,-buoy_prod)
plot(time,shear_prod - epsilon - buoy_prod,'--k','LineWidth',2)
% legend("d(tke)/dt","shear prod","- epsilon","- buoy prod","summed",'FontSize',20)
xlabel('time')
