% Script to check the implementation of div[psi]

addpath(genpath("../../../"));
clear; clc; close all;

%N = 100;
%for N = 50:100:1000

%% SOLUTION DOMAIN
% Cube size and position
L = 2*[1 1 1];
T = [5 5 3];

% Cube domain
%bndmesh = bndmeshCubeTranslated(N,L,T);

% Spherical domain
% bndmesh = mshSphere(N,1);
% bndmesh = bndmesh.translate(T);

%mesh = mshCube(N,L);
%mesh = mesh.translate(T);
%mesh = mesh.sub(1);
%bndmesh = mesh.bnd;

% Mesh size
%hvals(i) = sqrt(mean(bndmesh.ndv,1));

bndmesh = meshSymTetra();

Gamma = dom(bndmesh,12);

RWG = fem(bndmesh,'RWG');

KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;

N_gypsi = single_layer(Gamma,RWG.div,RWG.div);

Nelt = bndmesh.nelt;

[ii,jj] = meshgrid(1:Nelt,1:Nelt);

[Vel,DVel] = getTransVelDVel([1 0 0]);

N_SS = panel_assembly_shape_derivative(bndmesh,KV,RWG.div,RWG.div,ii(:),jj(:),Vel,DVel);

norm(N_gypsi-N_SS)/norm(N_SS)

P0 = fem(bndmesh,'P0');
SL = single_layer(Gamma,P0,P0);

delta = bndmesh.ndv;
delta = delta(1);

% 1,1 entry of N_gypsi in terms of SL
(2*SL(1,1)-2*SL(1,2))/delta^2

% For tetrahedron

%end