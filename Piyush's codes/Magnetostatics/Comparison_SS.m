% Comparing Galerkin Matrices with SS

addpath(genpath("../../"));
clear; clc; close all;

% Cube domain
bndmesh = bndmeshCubeTranslated(20,2*[1 1 1],[5 5 3]);

% Spherical domain
%bndmesh = mshSphere(N,1);
%bndmesh = bndmesh.translate([5 5 3]);

% Mesh size
%hvals(i) = sqrt(mean(bndmesh.ndv,1));

% BEM spaces
NED = fem(bndmesh,'NED'); % curl conforming -> Dirichlet trace
P1 = fem(bndmesh,'P1');
DIV0 = nxgrad(P1); % div conforming, div0 -> Neumann trace
DIV = fem(bndmesh,'RWG');

Gamma = dom(bndmesh,3);

Nelt = bndmesh.nelt;

