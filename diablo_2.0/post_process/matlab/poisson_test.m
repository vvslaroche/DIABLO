% test script for exploring FFT solve of Poisson equation
clc,clear,close all

NY = 11;
NZ = 9;
LY = 50;
LZ = 30;
% GY = zeros(1,NY);
% DY = (11:15) /10;
% 
% for j=1:NY
%     GY(j) = (j-1)*50/(NY-1) ;
%     GY(j) = (j-1)*50/(NY-1) ;
% end

GY = linspace(0,LY,NY);
GZ = linspace(0,LZ,NZ);


DY = (diff(GY(1:end-1))+diff(GY(2:end)))/2;
DY = [DY(1) DY DY(end)];
DZ = (diff(GZ(1:end-1))+diff(GZ(2:end)))/2;
DZ = [DZ(1) DZ DZ(end)];

% U = zeros(NY,NZ);
% VEC = U(:);
MATL = zeros(NY,NZ);
MATLL = zeros(NY,NZ);
MATD = zeros(NY,NZ);
MATU = zeros(NY,NZ);
MATUU = zeros(NY,NZ);
Mv1 = zeros(NY*NZ,NY*NZ);
Mv2 = zeros(NY*NZ,NY*NZ);

for j=2:NY-1
    for k=2:NZ-1
%         Mv2(j*k,k) = 
        MATL(j,k) = -1/(DY(j-1))^2;
        MATD(j,k) = 1/(DY(j))^2 + 1/(DY(j))^2  +  1/(DZ(k))^2 + 1/(DZ(k))^2;
        MATU(j,k) = -1/(DY(j+1))^2;
        MATLL(j,k) = -1/(DZ(k-1))^2;
        MATUU(j,k) = -1/(DZ(k+1))^2;
    end
end

for i1=1:NY*NZ
    for i2=1:NY*NZ
%         Mv1(i1,i2) = 
    end
end

for j=2:NY-1
    for k=2:NZ-1
        disp(['element j,k=' num2str(j) ',' num2str(k) ' is number ' num2str((j-1)*NZ+k)])
%         disp((j-1)*NZ+k)
        Mv1((j-1)*NZ+k,(j-1)*NZ+k) = MATD(j,k);
        Mv1((j-1)*NZ+k,(j-1)*NZ+k+1) = MATU(j,k);
        Mv1((j-1)*NZ+k,(j-1)*NZ+k-1) = MATL(j,k);
        Mv1((j-1)*NZ+k,(j-1)*NZ+k-1) = MATL(j,k);

%         Mv2(jk2i(j,k,NZ),jk2i(j,k,NZ)) = Mv2(jk2i(j,k,NZ),jk2i(j,k,NZ))+ 1/(DY(j))^2 + 1/(DY(j))^2  +  1/(DZ(k))^2 + 1/(DZ(k))^2;
%         Mv2(jk2i(j,k,NZ),jk2i(j-1,k,NZ)) = Mv2(jk2i(j,k,NZ),jk2i(j-1,k,NZ)) -1/(DY(j-1))^2;
%         Mv2(jk2i(j,k,NZ),jk2i(j+1,k,NZ)) = Mv2(jk2i(j,k,NZ),jk2i(j+1,k,NZ)) -1/(DY(j+1))^2;
%         Mv2(jk2i(j,k,NZ),jk2i(j,k-1,NZ)) = Mv2(jk2i(j,k,NZ),jk2i(j,k-1,NZ)) -1/(DZ(k-1))^2;
%         Mv2(jk2i(j,k,NZ),jk2i(j,k+1,NZ)) = Mv2(jk2i(j,k,NZ),jk2i(j,k+1,NZ)) -1/(DZ(k+1))^2;
    end
end

for j=1:NY
    for k=1:NZ
        disp(['element j,k=' num2str(j) ',' num2str(k) ' is number ' num2str((j-1)*NZ+k)])

        Mv2(jk2i(j,k,NZ),jk2i(j,k,NZ)) = Mv2(jk2i(j,k,NZ),jk2i(j,k,NZ))+ 1/(DY(j))^2 + 1/(DY(j))^2  +  1/(DZ(k))^2 + 1/(DZ(k))^2;
        if j>1% && j<5
            Mv2(jk2i(j,k,NZ),jk2i(j-1,k,NZ)) = Mv2(jk2i(j,k,NZ),jk2i(j-1,k,NZ)) -1/(DY(j-1))^2;
        end
        if j<NY% && j<5
            Mv2(jk2i(j,k,NZ),jk2i(j+1,k,NZ)) = Mv2(jk2i(j,k,NZ),jk2i(j+1,k,NZ)) -1/(DY(j+1))^2;
        end
        if k>1% && k>5
            Mv2(jk2i(j,k,NZ),jk2i(j,k-1,NZ)) = Mv2(jk2i(j,k,NZ),jk2i(j,k-1,NZ)) -1/(DZ(k-1))^2;
        end
        if k<NZ% && k>5
            Mv2(jk2i(j,k,NZ),jk2i(j,k+1,NZ)) = Mv2(jk2i(j,k,NZ),jk2i(j,k+1,NZ)) -1/(DZ(k+1))^2;
        end

        MATD(j,k) = 1/(DY(j))^2 + 1/(DY(j))^2  +  1/(DZ(k))^2 + 1/(DZ(k))^2;
        if j>1
            MATL(j,k) = -1/(DY(j-1))^2;
        end
        if j<NY
            MATU(j,k) = -1/(DY(j+1))^2;
        end
        if k>1
            MATLL(j,k) = -1/(DZ(k-1))^2;
        end
        if k<NZ
            MATUU(j,k) = -1/(DZ(k+1))^2;
        end
%         MATL(j,k) = -1/(DY(j-1))^2;
%         MATU(j,k) = -1/(DY(j+1))^2;
%         MATLL(j,k) = -1/(DZ(k-1))^2;
%         MATUU(j,k) = -1/(DZ(k+1))^2;
    end
end
% for k=1:NZ
%     for j=1:NY
%         j*k
%         Mv1(j*k,j*k) = MATD(j,k);
%     end
% end
% jk2i(2,2,NZ)


bandwidth = 2*NZ+1;
A = Mv2;
U = zeros(size(Mv2));
L = zeros(size(Mv2));

[truL,truU] = lu(A);

ops = 0;
% for i1=1:NY*NZ
for i1=1:NY*NZ
    L(i1,i1) = 1;

    lowerbound = 1;
    if i1>NZ+1
        lowerbound = i1-NZ;
    end
    upperbound = NY*NZ;
    if i1<NY*NZ-NZ
        upperbound = i1+NZ;
    end

%     for i2=1:NY*NZ
    for i2=lowerbound:upperbound
%         if i1>1
%             L(i1,i2) = A(i1,i2)/U(i2,i2);
%         end
        if i1>i2
            L(i1,i2) = A(i1,i2)/U(i2,i2);
            for i3=1:i2-1
                L(i1,i2) = L(i1,i2) - U(i3,i2)*L(i1,i3)/U(i2,i2);
                ops = ops+1;
            end
        end

        if i2>=i1
            U(i1,i2) = A(i1,i2);
            for i3=1:i1-1
                U(i1,i2) = U(i1,i2) - U(i3,i2)*L(i1,i3);
                ops = ops+1;
            end
        end

    end
%     if i1==25; break; end
end


b = rand(NY*NZ,1);
y = zeros(size(b));

truy = L\b;

for i1=1:NY*NZ
    y(i1) = b(i1);
    lowerbound = 1;
    if i1>NZ+1
        lowerbound = i1-NZ;
    end
    if i1>1
        for i2=lowerbound:i1-1
            y(i1) = y(i1) - L(i1,i2)*y(i2);
        end
    end
    y(i1) = y(i1)/L(i1,i1);
end



x = zeros(size(b));

trux = U\y;
% trux2 = A\b;

for i1=NY*NZ:-1:1
    x(i1) = y(i1);
    upperbound = NY*NZ;
    if i1<NY*NZ-NZ
        upperbound = i1+NZ;
    end
    if i1<NY*NZ
        for i2=i1+1:upperbound
            x(i1) = x(i1) - U(i1,i2)*x(i2);
        end
    end
    x(i1) = x(i1)/U(i1,i1);
end


function i=jk2i(j,k,NZ)
i=(j-1)*NZ+k;
end