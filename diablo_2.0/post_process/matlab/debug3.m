% code that test different combinations of grid sizes to see the size of
% global array variables

clc,clear

NX = 32;
NY = 65;
NYP = 33;
NZ = 49;
NZP = 25;

NPROCZ = 2;
NXP = NX/(2*NPROCZ);


disp(' U size: ')
disp((NX+2)*(NZP+2)*(NYP+2))
disp('CU size: ')
disp(2*(NXP+1)*(NZ+2)*(NYP+2))


%% channel

NX = 16:300;
NZ = 16:256;
NPROCZ = 16;
U = zeros(length(NX),length(NZ));
CU = zeros(length(NX),length(NZ));

NZP = NZ/NPROCZ;
NXP = NX/(2*NPROCZ);

for i=1:length(NX)
    for j=1:length(NZ)
        U(i,j) = (NX(i)+2)*(NZP(j)+2);
        CU(i,j) = 2*(NXP(i)+1)*(NZ(j)+2);
    end
end
r = (U>CU)*1.0;

figure
pcolor(NX,NZ,r')
colorbar
shading flat
xlabel('NX'), ylabel('NZ')
