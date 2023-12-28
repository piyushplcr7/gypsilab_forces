% New transmission problem script

addpath(genpath("../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 4;
mu0 = 2;
vals = 5:10;
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
forces_bem = forces_mst;
torques_mst = forces_mst;
torques_bem = forces_mst;
hvals = 0*vals;

for i = 1:Nvals
    N = 2^vals(i);
    disp(N);
    %% SOLUTION DOMAIN
    % Cube size and position
    L = 2*[1 1 1];
    T = 0*[1 0.5 2];
    
    % Bounding box
    bndmesh_e = mshSphere(N,5);

    % Spherical domain
%     bndmesh_i = mshSphere(N,1);
%     bndmesh_i = bndmesh_i.translate(T);
    
    mesh = mshCube(N,L);
    mesh = mesh.translate(T);
    %mesh = mesh.sub(1);
    bndmesh_i = mesh.bnd;
    
    % Mesh size
    hvals(i) = sqrt(mean(bndmesh_i.ndv,1));

    % Visualize with outer mesh translucent
    %customPlot(bndmesh_i,bndmesh_e);
    
    Gamma_i = dom(bndmesh_i,3);
    normals_i = Gamma_i.qudNrm;
    
    %% Solving the transmission problem
    % These are traces from the exterior
    B_0 = [1,1,1];
    [Psi_i,g_i,Psi_e] = solveTPConstantFld(bndmesh_i,bndmesh_e,mu,mu0,B_0);
    
    %% MST Based force and torque on bndmesh_i
    % Force computation
    NED_i = fem(bndmesh_i,'NED'); 
    P1_i = fem(bndmesh_i,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0_i = nxgrad(P1_i); 
    RWG_i = fem(bndmesh_i,'RWG');

    % Bn = curlA.n = curlTg
    Bn = reconstruct(g_i,Gamma_i,NED_i.curl);
    % Htangential = nx(Hxn) = mu_0^-1 nxPsi
    % Psi_i is restriction Neumann trace from \partial\Omega_e to \Gamma_i

    Psivals_i = reconstruct(Psi_i,Gamma_i,DIV0_i);
    Ht = 1/mu0 * cross(normals_i,Psivals_i,2);
    Ht = vecnorm(Ht,2,2);

    forces_mst(i,:) = ForceMstTP(Gamma_i,Bn,Ht,mu0,mu)

    % Torque computation
    Xcg = [4 0 0];
    torques_mst(i,:) = TorqueMstTP(Gamma_i,Bn,Ht,mu0,mu,Xcg)
end
