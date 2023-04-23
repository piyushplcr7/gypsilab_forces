% New transmission problem script

addpath(genpath("../../../"));
clear; clc; close all;

% (mui+mue)/(mui-mue)
mu = 3;
mu0 = 2;
vals = 50:100:1000;
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
forces_bem = forces_mst;
forces_bem_1 = forces_mst;

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
    
    %% Solving the transmission problem
    % These are traces from the exterior
    % Psi lies in the space nxgradP1 and g lies in the space NED
    [Psi,g] = solveTransmissionProblem(bndmesh,J,omega_src,mu,mu0);
    % Interior traces
    Psi_in = mu/mu0 * Psi;
    g_in = g;
    
    %% Computing the MST based force and torque
    
    % Force computation
    NED = fem(bndmesh,'NED'); 
    P1 = fem(bndmesh,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0 = nxgrad(P1); 
    RWG = fem(bndmesh,'RWG');
    
    % Bn = curlA.n = curlTg
    Bn = reconstruct(g,Gamma,NED.curl);
    % Ht = nx(Hxn) = mu_e^-1 nxPsi
    Psivals = reconstruct(Psi,Gamma,DIV0);
    Ht = mu0^(-1) * cross(normals,Psivals,2);
    Ht = vecnorm(Ht,2,2);
    
    forces_mst(i,:) = ForceMstTP(Gamma,Bn,Ht,mu0,mu)

    % Torque computation

    %% Computing forces and torques using BEM shape derivative

    % Projecting Psi to RWG
    Psi_RWG = proj(Psivals,Gamma,RWG);

    [Vel1,DVel1] = getTransVelDVel([1 0 0]);
    [Vel2,DVel2] = getTransVelDVel([0 1 0]);
    [Vel3,DVel3] = getTransVelDVel([0 0 1]);

%     f1 = -1/mue * SuperConductorShapeDerivativeT3(bndmesh,Psi_RWG,Vel1,omega_src,J)...
%         + 1/mue * ShapeDerivativel2BEMTP(bndmesh,g,Vel1,omega_src,J);
%     f2 = -1/mue * SuperConductorShapeDerivativeT3(bndmesh,Psi_RWG,Vel2,omega_src,J)...
%         + 1/mue * ShapeDerivativel2BEMTP(bndmesh,g,Vel2,omega_src,J);
%     f3 = -1/mue * SuperConductorShapeDerivativeT3(bndmesh,Psi_RWG,Vel3,omega_src,J)...
%         + 1/mue * ShapeDerivativel2BEMTP(bndmesh,g,Vel3,omega_src,J);
% 
%     forces_bem(i,:) = [f1 f2 f3]

    nf1 = ForceSdBemTP(bndmesh,Psi,g,J,omega_src,Vel1);
    nf2 = ForceSdBemTP(bndmesh,Psi,g,J,omega_src,Vel2);
    nf3 = ForceSdBemTP(bndmesh,Psi,g,J,omega_src,Vel3);

    forces_bem_1(i,:) = [nf1 nf2 nf3]


end
