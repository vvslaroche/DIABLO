% calculates divergence using 3D flow field for duct flow

clear,clc,close all

datapath = '/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-data/';
dirname = 'duct-KH-NY256-NU1.0e-4';

SAVE_FLOW_INT = 8048;

NX = 192;
NY = 257;
NZ = 129;

LX = 14;
gx = linspace(0,LX,NX)';
gxf = (gx(2:end)+gx(1:end-1))/2;

% LZ = 10;
% gx = linspace(0,LZ,NZ)';
% gxf = (gx(2:end)+gx(1:end-1))/2;

filename = [datapath dirname '/grid.h5'];
gy = h5read(filename,'/grids/y');
gz = h5read(filename,'/grids/z');
gyf = (gy(1:end-1)+gy(2:end))/2;
gzf = (gz(1:end-1)+gz(2:end))/2;
LY = gyf(end)-gyf(1);
LZ = gzf(end)-gzf(1);


for i=1:200
    filename=[datapath dirname '/output/out' num2str(SAVE_FLOW_INT*i,'%05.0f') '.h5'];
    disp(['Reading out' num2str(SAVE_FLOW_INT*i,'%05.0f') '.h5'])
    try
        U = h5read(filename,'/Timestep/U');
    catch
        disp(['out' num2str(SAVE_FLOW_INT*i,'%05.0f') '.h5 does not exist. Exiting...'])
        break
    end
    V = h5read(filename,'/Timestep/V');
    W = h5read(filename,'/Timestep/W');
    dudx = (U(2:end,:,:) - U(1:end-1,:,:)) ./ reshape((gx(2:end) -gx(1:end-1)) ,[length(gx)-1,1,1]);
%     dudx = 0;
    dvdy = (V(:,2:end,:) - V(:,1:end-1,:)) ./ reshape((gyf(2:end)-gyf(1:end-1)),[1,length(gyf)-1,1]);
    dwdz = (W(:,:,2:end) - W(:,:,1:end-1)) ./ reshape((gzf(2:end)-gzf(1:end-1)),[1,1,length(gzf)-1]);
%     figure, imagesc(squeeze(dvdy(1,:,:))), colorbar, title('dvdy test')
%     figure, imagesc(squeeze(dwdz(1,:,:))), colorbar, title('dwdz test')
    
    % interpolate off-derivative directions onto the trimmed grid
    dvdyN = interp1(gzf,permute(dvdy,[3,2,1]),gz(2:end-1));
    dvdyN = permute(dvdyN,[3,2,1]);
    dvdyN = interp1(gx,dvdyN,gxf);

    dwdzN = interp1(gyf,permute(dwdz,[2,1,3]),gy(2:end-1));
    dwdzN = permute(dwdzN,[2,1,3]);
    dwdzN = interp1(gx,dwdzN,gxf);

    dudxN = interp1(gyf,permute(dudx,[2,1,3]),gy(2:end-1));
    dudxN = permute(dudxN,[2,1,3]);
    dudxN = interp1(gzf,permute(dudxN,[3,2,1]),gz(2:end-1));
    dudxN = permute(dudxN,[3,2,1]);
    
    div = dvdyN + dwdzN + dudxN;
%     figure, imagesc(squeeze(dvdyN(1,:,:))), colorbar, title('dvdy')
%     figure, imagesc(squeeze(dwdzN(1,:,:))), colorbar, title('dwdz')
%     figure, imagesc(squeeze(div(1,:,:))), colorbar, title('div')
    figure, imagesc(squeeze(div(1,2:end-1,2:end-1))), colorbar, title('div')
%     clim([-1e-1 1e-1])
    
    % div = sum(dudx,'all') + sum(dvdy,'all') + sum(dwdz,'all');
    % div = trapz(gzf,trapz(gy(2:end-1),abs(dvdy))) ...
    %     + trapz(gyf,trapz(gz(2:end-1),abs(dwdz)'));

    divavg = trapz(gz(2:end-1), ...
             trapz(gy(2:end-1), ...
             trapz(gxf,div))) /LX/LY/LZ;
%     divavg = trapz(gz(3:end-2), ...
%              trapz(gy(3:end-2), ...
%              trapz(gxf,div(:,2:end-1,2:end-1)))) /LX/LY/LZ;
    disp(divavg)
end

