% Script testing CUDA
addpath(genpath("../../../"));

clear;clc;

format long;

%bndmesh = mshSphere(2,1);
mesh = mshCube(2,[1 1 1]);
% mesh = mesh.sub(1);
bndmesh = mesh.bnd;

% bndmesh = meshSymTetra;

P0 = fem(bndmesh,'P0');
P1 = fem(bndmesh,'P1');
RWG = fem(bndmesh,'RWG');

KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
KK = @(x,y,z) -z./vecnorm(z,2,2).^3/4./pi;
Nelt = bndmesh.nelt;
[ii,jj] = meshgrid(1:Nelt,1:Nelt);

tic;
% V = panel_assembly(bndmesh,KV,P0,P0,ii(:),jj(:));
% V = panel_assembly(bndmesh,KK,RWG,RWG,ii(:),jj(:));
% V = panel_assembly(bndmesh,KK,ntimes(P1),P0,ii(:),jj(:));
V = panel_assembly(bndmesh,KV,P1,P1,ii(:),jj(:));


elapsed_time = toc;
fprintf('matlab took %.4f seconds.\n', elapsed_time);


tic;
V_mex = getSSMatrixMex(bndmesh);
elapsed_time = toc;
fprintf('C++ took %.4f seconds.\n', elapsed_time);

norm(V-V_mex)
norm(V-V_mex)/norm(V)