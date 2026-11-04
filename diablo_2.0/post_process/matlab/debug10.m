% calculates divergence using 3D flow field for CHANNEL flow

clear,clc,close all

% datapath = '/Users/vincent/Desktop/';
% dirname = 'spectra_test';
datapath = '/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-data/';
dirname = 'duct-chan-KH-hpc-highres';

% SAVE_FLOW_INT = 500;
SAVE_FLOW_INT = 200;

% NX = 768;
% NY = 1153;
% NZ = 256;
NX = 96;
NY = 129;
NZ = 64;

LX = 13.98;
gx = linspace(0,LX,NX)';
gxf = (gx(2:end)+gx(1:end-1))/2;
% LZ = 5;
LZ = 10;
gz = linspace(0,LZ,NZ)';
gzf = (gz(2:end)+gz(1:end-1))/2;

filename = [datapath dirname '/grid.h5'];
gy = h5read(filename,'/grids/y');
gyf = (gy(1:end-1)+gy(2:end))/2;
LY = gyf(end)-gyf(1);


for i=1:200
    filename=[datapath dirname '/output/out' num2str(SAVE_FLOW_INT*i,'%05.0f') '.h5'];
    disp(['Reading out' num2str(SAVE_FLOW_INT*i,'%05.0f') '.h5'])
    try
        U = h5read(filename,'/Timestep/U');
    catch
        disp(['out' num2str(SAVE_FLOW_INT*i,'%05.0f') '.h5 does not exist. Exiting...'])
        break
    end

    dudx = (U(2:end,:,:) - U(1:end-1,:,:)) ./ reshape((gx(2:end) -gx(1:end-1)) ,[length(gx)-1,1,1]);
    clear U

    dudx = interp1(gyf,permute(dudx,[2,1,3]),gy(2:end-1));
    dudx = permute(dudx,[2,1,3]);
    dudx = interp1(gz,permute(dudx,[3,2,1]),gzf);
    dudx = permute(dudx,[3,2,1]);

    div = dudx;
    clear dudx


    V = h5read(filename,'/Timestep/V');
    dvdy = (V(:,2:end,:) - V(:,1:end-1,:)) ./ reshape((gyf(2:end)-gyf(1:end-1)),[1,length(gyf)-1,1]);
    clear V

    dvdy = interp1(gz,permute(dvdy,[3,2,1]),gzf);
    dvdy = permute(dvdy,[3,2,1]);
    dvdy = interp1(gx,dvdy,gxf);

    div = div + dvdy;
    clear dvdy


    W = h5read(filename,'/Timestep/W');
    dwdz = (W(:,:,2:end) - W(:,:,1:end-1)) ./ reshape((gz(2:end) -gz(1:end-1)) ,[1,1,length(gz)-1]);
    clear W

    dwdz = interp1(gyf,permute(dwdz,[2,1,3]),gy(2:end-1));
    dwdz = permute(dwdz,[2,1,3]);
    dwdz = interp1(gx,dwdz,gxf);

    div = div + dwdz;
    clear dwdz
%     figure, imagesc(squeeze(dvdy(1,:,:))), colorbar, title('dvdy test')
%     figure, imagesc(squeeze(dwdz(1,:,:))), colorbar, title('dwdz test')
    
    % interpolate off-derivative directions onto the trimmed grid

% 
%     dwdz = interp1(gyf,permute(dwdz,[2,1,3]),gy(2:end-1));
%     dwdz = permute(dwdz,[2,1,3]);
%     dwdz = interp1(gx,dwdz,gxf);
% 
%     dudx = interp1(gyf,permute(dudx,[2,1,3]),gy(2:end-1));
%     dudx = permute(dudx,[2,1,3]);
%     dudx = interp1(gz,permute(dudx,[3,2,1]),gzf);
%     dudx = permute(dudx,[3,2,1]);
    
%     div = dvdy + dwdz + dudx;
%     figure, imagesc(squeeze(dvdyN(1,:,:))), colorbar, title('dvdy')
%     figure, imagesc(squeeze(dwdzN(1,:,:))), colorbar, title('dwdz')
    figure, imagesc(squeeze(div(:,:,1))), colorbar, title('div')
    clim([-1e-1 1e-1])
    
    % div = sum(dudx,'all') + sum(dvdy,'all') + sum(dwdz,'all');
    % div = trapz(gzf,trapz(gy(2:end-1),abs(dvdy))) ...
    %     + trapz(gyf,trapz(gz(2:end-1),abs(dwdz)'));
    divavg = trapz(gzf, ...
             trapz(gy(2:end-1), ...
             trapz(gxf,div))) /LX/LY/LZ;
%     divsum = trapz(gz(2:end-1),trapz(gy(2:end-1),abs(squeeze(div(1,:,:)))));
    disp(divavg)
end

