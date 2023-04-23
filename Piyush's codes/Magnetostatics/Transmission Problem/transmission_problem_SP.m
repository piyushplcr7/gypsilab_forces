% New transmission problem script

addpath(genpath("../../../"));
clear; clc; close all;

% (mui+mue)/(mui-mue)
mu = 3;
mu0 = 2;
vals = 50:100:1000;
Nvals = size(vals,2);
forces_vol = zeros(Nvals,3);
forces_bem = forces_vol;
forces_bem_1 = forces_vol;

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
    
    %mesh = mshCube(N,L);
    %mesh = mesh.translate(T);
    %mesh = mesh.sub(1);
    %bndmesh = mesh.bnd;
    
    % Mesh size
    %hvals(i) = sqrt(mean(bndmesh.ndv,1));
    
    Gamma = dom(bndmesh,3);
    normals = Gamma.qudNrm;
    
    %% Source
    N_src = N;
    R0 = 2;
    r0 = .5;
    [J,mesh_src] = get_torus_source(N_src,R0,r0);
    omega_src = dom(mesh_src,3);
    
    %% Solving the transmission problem using scalar potential formulation
    % These are traces from the exterior
    % Psi lies in the space nxgradP1 and g lies in the space NED
    [Psi,g] = solveTpScalPotBIE(bndmesh,mu,mu0,J,omega_src);
    
    %% Computing the MST based force and torque
    [Vel1,DVel1] = getTransVelDVel([1 0 0]);
    [Vel2,DVel2] = getTransVelDVel([0 1 0]);
    [Vel3,DVel3] = getTransVelDVel([0 0 1]);

    f1 = shapeDer_ScalPot_Vol_TP(bndmesh,mu,mu0,Psi,g,J,omega_src,Vel1,DVel1);
    f2 = shapeDer_ScalPot_Vol_TP(bndmesh,mu,mu0,Psi,g,J,omega_src,Vel2,DVel2);
    f3 = shapeDer_ScalPot_Vol_TP(bndmesh,mu,mu0,Psi,g,J,omega_src,Vel3,DVel3);

    forces_vol(i,:) = [f1 f2 f3]


end
