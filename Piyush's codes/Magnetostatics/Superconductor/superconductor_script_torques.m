% Superconductor script

addpath(genpath("../../../"));
clear; clc; close all;
format long;
mui = 1;
mue = 1;
%N = 50;

for N=50:100:1000
disp(N)
%for N=50
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

%plot_field(TnA,bndmesh,J,omega_src);

%% Visualizing the MST force density 
% figure;
% %plot(bndmesh);
% %hold on;
% % Coefficients for zero Dirichlet Trace
% TdA = TnA * 0;
% Gamma = dom(bndmesh,3);
% VisualizeForceDensity(TdA,TnA,Gamma);

%% Computing forces
% Dom object
Gamma = dom(bndmesh,3);

% Coefficients for zero Dirichlet Trace
TdA = TnA * 0;

Xcg = [4 0 0];
torque_mst = MstTorqueFromA(TdA,TnA,Gamma,Xcg)

% Shape Derivative computation of torque
% Getting Rotational Vels and DVels
% [Velxr,DVelxr] = getRotVelDVel([1 0 0],Xcg);
% [Velyr,DVelyr] = getRotVelDVel([0 1 0],Xcg);
% [Velzr,DVelzr] = getRotVelDVel([0 0 1],Xcg);
% 
% sd_e1 = SuperConductorShapeDerivative(bndmesh,TnA,Velxr,DVelxr,omega_src,J);
% sd_e2 = SuperConductorShapeDerivative(bndmesh,TnA,Velyr,DVelyr,omega_src,J);
% sd_e3 = SuperConductorShapeDerivative(bndmesh,TnA,Velzr,DVelzr,omega_src,J);
% 
% torque_sd = [sd_e1 sd_e2 sd_e3]

save("SC_VP_Torques.mat","forces_mst","forces_sd","torques_mst","torques_bem","hvals");

end