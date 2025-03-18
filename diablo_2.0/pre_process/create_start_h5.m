clc, clear

filename = 'start.h5';
NX = 128;
NY = 193; % total NY
NZ = 32;
N_TH = 7;
KICK = 0.1;
t_start = 0.0; % Start time
R = 3;
LX = 15.22;  % 13.98, 15.22
LZ = 5;

gy = h5read('grid.h5','/grids/y')';
gyf = (gy(1:end-1)+gy(2:end))/2;
gx = reshape(linspace(0,LX,NX+1),[NX+1 1 1]);
gz = reshape(linspace(0,LZ,NZ+1),[1 1 NZ+1]);



% Create flow fields
V = zeros(NX, NY, NZ);
W = zeros(NX, NY, NZ);
TH = zeros(NX, NY, NZ, N_TH);

U = tanh(gyf).*ones(NX,1,1).*ones(1,1,NZ);
TH(:,:,:,1) = tanh(R*gyf)  .*ones(NX,1,1).*ones(1,1,NZ);
% TH(:,:,:,2) = (tanh(gyf-1)+1)/2  .*ones(NX,1,1).*ones(1,1,NZ);
% TH(:,:,:,3) = (-tanh(gyf+1)+1)/2  .*ones(NX,1,1).*ones(1,1,NZ);


load('/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-3D/Unstable_Phytoplankton/SS_h07.1.mat');
% load('/Users/vincent/Documents/_School/Berkeley/Research/DIABLO-3D/Unstable_Phytoplankton/SS_h10.0.mat');
newP = interp1(y+15,P,gyf);
newN = interp1(y+15,N,gyf);
TH(:,:,:,2) = newP .*ones(NX,1,1).*ones(1,1,NZ);
TH(:,:,:,3) = newN .*ones(NX,1,1).*ones(1,1,NZ);
TH(:,:,:,4) = TH(:,:,:,2);
TH(:,:,:,5) = TH(:,:,:,3);
TH(:,:,:,6) = TH(:,:,:,2);
TH(:,:,:,7) = TH(:,:,:,3);



% add random noise
s = RandStream("dsfmt19937");
U = U+KICK*(rand(s,size(U))-0.5);
V = V+KICK*(rand(s,size(V))-0.5);
W = W+KICK*(rand(s,size(W))-0.5);


% add primary mode perturbation
U = U -0.01*2*cos(2*pi*gx(1:NX)/LX).*tanh(gyf)./cosh(gyf) .*ones(1,1,NZ);
V = V +0.01*2*2*pi/LX*sin(2*pi*gx(1:NX)/LX)./cosh(gyf) .*ones(1,1,NZ);


% write to file
h5create(filename, '/Timestep/U', size(U));
h5write(filename, '/Timestep/U', U);
h5create(filename, '/Timestep/V', size(V));
h5write(filename, '/Timestep/V', V);
h5create(filename, '/Timestep/W', size(W));
h5write(filename, '/Timestep/W', W);
for n=1:N_TH
    h5create(filename, ['/Timestep/TH' num2str(n,'%2.2d')], size(TH(:,:,:,n)));
    h5write(filename, ['/Timestep/TH' num2str(n,'%2.2d')], TH(:,:,:,n));
end

h5writeatt(filename, '/Timestep/','Time', t_start);
h5writeatt(filename, '/', 'Resolution', [NX; NY; NZ]);

