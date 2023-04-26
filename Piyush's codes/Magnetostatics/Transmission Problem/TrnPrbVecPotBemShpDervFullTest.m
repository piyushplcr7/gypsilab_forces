% Full BEM shape derivative Transmission Problem Vect

addpath(genpath("../../../"));
clear; clc; close all;

% (mui+mue)/(mui-mue)
mu = 4;
mu0 = 2;
vals = 50:100:1000;
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
forces_bem = forces_mst;
forces_bem_simplified = forces_mst;
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

    %% Comparing forces from full Shape Derivative and simplified Shape Derivative
    NED = fem(bndmesh,'NED'); 
    P1 = fem(bndmesh,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0 = nxgrad(P1); 
    RWG = fem(bndmesh,'RWG');
    
    % Projecting Psi to RWG
    Psivals = reconstruct(Psi,Gamma,DIV0);
    Psi_RWG = proj(Psivals,Gamma,RWG);

    [Vel1,DVel1] = getTransVelDVel([1 0 0]);
    [Vel2,DVel2] = getTransVelDVel([0 1 0]);
    [Vel3,DVel3] = getTransVelDVel([0 0 1]);

    nf1 = ForceSdBemTP(bndmesh,Psi,g,J,omega_src,Vel1);
    nf2 = ForceSdBemTP(bndmesh,Psi,g,J,omega_src,Vel2);
    nf3 = ForceSdBemTP(bndmesh,Psi,g,J,omega_src,Vel3);

    forces_bem_simplified(i,:) = [nf1 nf2 nf3]

    sd1 = SdBemTPVP(bndmesh,Psi_RWG,g,J,omega_src,Vel1,DVel1,mu0,mu);
    sd2 = SdBemTPVP(bndmesh,Psi_RWG,g,J,omega_src,Vel2,DVel2,mu0,mu);
    sd3 = SdBemTPVP(bndmesh,Psi_RWG,g,J,omega_src,Vel3,DVel3,mu0,mu);

    forces_bem(i,:) = [sd1 sd2 sd3]


end
or Potential test