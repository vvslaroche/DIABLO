clear

NY=11;
NZ=5;

U1 = ones(NY,NZ);
U2 = ones(NY,NZ);
U3 = ones(NY,NZ);
VEC = zeros(NY,NZ);
MATL = zeros(NY,NZ);
MATD = zeros(NY,NZ);
MATU = zeros(NY,NZ);


for j=1:NY
    for k=1:NZ
        VEC(j,k) = U1(j,k);
    end
end


%%
loc=zeros(NY,NZ);
for j=1:NY
    for k=1:NZ
        loc(j,k) = (k-1)*NY+j;

    end
end
VEC = loc(:);
MAT = zeros(NY*NZ);

DY=5;
DZ=3;

% for j=1:NY
%     for k=1:NZ
%         MAT(loc(j,k),loc(j,k)) = -2/DY -2/DZ;        % j, k 
%         if j~=NY
%             MAT(loc(j,k),loc(j+1,k)) = 1/DY;         % j+1, k
%         end
%         if j~=1
%             MAT(loc(j,k),loc(j-1,k)) = 1/DY;         % j-1, k
%         end
%         if k~=NZ
%             MAT(loc(j,k),loc(j,k+1)) = 1/DZ;         % j, k+1
%         end
%         if k~=1
%             MAT(loc(j,k),loc(j,k-1)) = 1/DZ;         % j, k-1
%         end
%     end
% end
for j=2:NY-1
    for k=2:NZ-1
        MAT(loc(j,k),loc(j,k)) = -2/DY -2/DZ;        % j, k
        MAT(loc(j,k),loc(j+1,k)) = 1/DY;         % j+1, k
        MAT(loc(j,k),loc(j-1,k)) = 1/DY;         % j-1, k
        MAT(loc(j,k),loc(j,k+1)) = 1/DZ;         % j, k+1
        MAT(loc(j,k),loc(j,k-1)) = 1/DZ;         % j, k-1
    end
end