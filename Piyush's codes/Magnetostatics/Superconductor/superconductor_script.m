% Superconductor script

addpath(genpath("../../../"));
clear; clc; close all;

mui = 1;
mue = 1;
%N = 50;

%for N=50:100:2000
for N=50
%% SOLUTION DOMAIN
% Cube size and position
% L = 2*[1 1 1];
 T = [5 5 3];

% Cube domain
% bndmesh = bndmeshCubeTranslated(N,L,T);

% Spherical domain
bndmesh = mshSphere(N,1);
bndmesh = bndmesh.translate(T);

%mesh = mshCube(N,L);
%mesh = mesh.translate(T);
%mesh = mesh.sub(1);
%bndmesh = mesh.bnd;

% Mesh size
%hvals(i) = sqrt(mean(bndmesh.ndv,1));

%% Current source
N_src = N;
R0 = 2;
r0 = .5;
[J,mesh_src] = get_torus_source(N_src,R0,r0);
omega_src = dom(mesh_src,3);

%% Solving the problem and obtaining the Neumann trace
[TnA0,TnA] = solve_superconductor(bndmesh,J,mesh_src);

%% Plotting the computed B field

plot_field(TnA,bndmesh,J,omega_src);

%% Computing forces
% Dom object
Gamma = dom(bndmesh,3);

% Coefficients for zero Dirichlet Trace
TdA = TnA * 0;

force_mst = MstForceFromA(TdA,TnA,Gamma)
end