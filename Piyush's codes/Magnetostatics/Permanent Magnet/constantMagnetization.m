% Script for constant Magnetization

addpath(genpath("../../../"));
clear; clc; close all;
format long;

mu = 1;
mu0 = 1;
vals = 50:100:1000;
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
    J = @(x,y,z) 1 * J_orig(x,y,z);

    % These are traces from the exterior
    % Psi lies in the space nxgradP1 and g lies in the space NED
    [TnAJ,TdAJ,TnAM,TdAM] = solveMagnetProblemSimplified(Gamma,J,omega_src,mu,mu0,M);

    % Projecting TnAM from nxgradP1 to RWG
    P1 = fem(bndmesh,'P1');
    RWG = fem(bndmesh,'RWG');
    PsivalsM = reconstruct(TnAM,Gamma,nxgrad(P1));
    TnAM_RWG = proj(PsivalsM,Gamma,RWG);

    PsivalsJ = reconstruct(TnAJ,Gamma,nxgrad(P1));
    TnAJ_RWG = proj(PsivalsJ,Gamma,RWG);

    %% Computing the B field for verification
%     plot_field_magnet(TdAM,TnAM_RWG,bndmesh,J,omega_src,mu0,interior);
%     figure;
%     plot_field_magnet(TdAJ,TnAJ_RWG,bndmesh,J,omega_src,mu0,interior);

    %% Computing forces and torques from MST

    % Force computation
    NED = fem(bndmesh,'NED'); 
    P1 = fem(bndmesh,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0 = nxgrad(P1); 
    RWG = fem(bndmesh,'RWG');
    
    % BJn = curlAJ.n = curlTgJ
    BJn = reconstruct(TdAJ,Gamma,NED.curl);
    % HJt = nx(HJxn) = mu0^-1 nxPsiJ
    PsiJvals = reconstruct(TnAJ,Gamma,DIV0);
    HJt = mu0^(-1) * cross(normals,PsiJvals,2);
    HJt = vecnorm(HJt,2,2);

    % Traces for the magnetic part
    BMn = reconstruct(TdAM,Gamma,NED.curl);
    % curlAM x n = BM x n = mu0 HM x n from outside
    PsiMvals = reconstruct(TnAM,Gamma,DIV0);
    % Tangential component of HM
    HMt = 1/mu0 * cross(normals,PsiMvals,2);
    %HMt = vecnorm(HMt,2,2);
    
    forces_mst(i,:) = -ForceMstMagnet(Gamma,BMn,HMt,mu0,mu,M)+ForceMstTP(Gamma,BJn,HJt,mu0,mu)

    % Torque computation
    Xcg = [4 0 0];
    torques_mst(i,:) = -TorqueMstMagnet(Gamma,BMn,HMt,mu0,mu,M,Xcg)+TorqueMstTP(Gamma,BJn,HJt,mu0,mu,Xcg)

end
