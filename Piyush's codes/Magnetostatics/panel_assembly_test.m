% Panel Assembly Test Script
addpath(genpath("../../"));
clear; clc; close all;
% Computing the Galerkin Matrix Corresponding to A operator

% Gypsi Computation

N = 50;
T = [5 5 3];

% Cube domain
% bndmesh = bndmeshCubeTranslated(N,L,T);

% Spherical domain
bndmesh = mshSphere(N,1);
bndmesh = bndmesh.translate(T);

% Div conforming space 
DIV = fem(bndmesh,'RWG');

% Dom objects for integration
Gamma = dom(bndmesh,3);

% Galerkin matrix for A
Amat_gypsi = single_layer(Gamma,DIV,DIV);

% Panel Assembly computation
% Defining the kernel for single layer BIO
% KV = @(x,y,z) 1/norm(z)/4./pi;
KV = @(x,y,z) sqrt(1./ sum(z.^2 ,2) ) /4./pi;

Nelt = bndmesh.nelt;

[ii,jj] = meshgrid(1:Nelt,1:Nelt);

Amat_panel_assembly = panel_assembly(bndmesh,KV,DIV,DIV,ii(:),jj(:));



