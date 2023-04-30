% Script for constant Magnetization

addpath(genpath("../../../"));
clear; clc; close all;
format long;

mu = 1;
mu0 = 1;
vals = 250;
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
forces_bem = forces_mst;
torques_mst = forces_mst;
torques_bem = forces_mst;

for i = 1:Nvals
    N = vals(i);
    disp(N);
    %% SOLUTION DOMAIN
    % Cube size and position
    L = 2*[1 1 1];
    T = [5 5 3];
    
    % Cube domain
    %bndmesh = bndmeshCubeTranslated(N,L,T);
    
    % Spherical domain
    bndmesh = mshSphere(N,1);
    bndmesh = bndmesh.translate(T);
    
%     mesh = mshCube(N,L);
%     mesh = mesh.translate(T);
    %mesh = mesh.sub(1);
%     bndmesh = mesh.bnd;
    
    % Mesh size
    %hvals(i) = sqrt(mean(bndmesh.ndv,1));
    
    Gamma = dom(bndmesh,3);
    normals = Gamma.qudNrm;

    % Function to determine the interior of the magnet
    interior = @(X) (vecnorm(X-T,2,2) < 1.05);
    
    %% Source
    N_src = N;
    R0 = 2;
    r0 = .5;
    [J_orig,mesh_src] = get_torus_source(N_src,R0,r0);
    omega_src = dom(mesh_src,3);
    
    %% Solving the Magnet problem
    % Constant Magnetization function
    M = @(X) ones(size(X,1),1) * [1 0 0];
    % Modifying J source
    %J = @(x,y,z) 0 * J_orig(x,y,z);

    % These are traces from the exterior
    % Psi lies in the space nxgradP1 and g lies in the space NED
    [TnA,TdA] = solveMagnetProblemSimplified(Gamma,J,omega_src,mu,mu0,M);

    % Projecting TnA from nxgradP1 to RWG
    P1 = fem(bndmesh,'P1');
    RWG = fem(bndmesh,'RWG');
    Psivals = reconstruct(TnA,Gamma,nxgrad(P1));

    TnA_RWG = proj(Psivals,Gamma,RWG);

    %% Computing the B field for verification
    plot_field_magnet(TdA,TnA_RWG,bndmesh,J,omega_src,mu0,interior);
end
