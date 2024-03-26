% Script to compare the operator evaluated using gypsilab and mex file

N = 512;
L = 2*[1 1 1];
T = [5 5 3];

mesh = mshCube(N,L);
mesh = mesh.translate(T);
bndmesh = mesh.bnd;

P1 = fem(bndmesh,'P1');
P0 = fem(bndmesh,'P0');

Gamma = dom(bndmesh,3);

V = single_layer(Gamma,P0,P0);
K = double_layer_laplace(Gamma,P0,P1);
W = single_layer(Gamma,P1.nxgrad,P1.nxgrad);

%% Calling the mex function

%% Assembling using SS in gypsi
 Nelt = bndmesh.nelt;
[ii,jj] = meshgrid(1:Nelt,1:Nelt);

SLKernel = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
DLKernel = @(x,y,z) -1/4/pi * z./vecnorm(z,2,2).^3;

VSS_gypsi = panel_assembly(bndmesh,SLKernel,P0,P0,ii,jj);
KSS_gypsi = panel_assembly(bndmesh,DLKernel,ntimes(P1),P0,ii,jj);
WSS_gypsi = panel_assembly(bndmesh,SLKernel,P1.nxgrad,P1.nxgrad,ii,jj);





