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

KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
Nelt = bndmesh.nelt;
[ii,jj] = meshgrid(1:Nelt,1:Nelt);

% tic;
V = panel_assembly(bndmesh,KV,P0,P0,ii(:),jj(:));

% elapsed_time = toc;
% fprintf('CPU took %.4f seconds.\n', elapsed_time);

V_mex = getSSMatrixMex(bndmesh);
