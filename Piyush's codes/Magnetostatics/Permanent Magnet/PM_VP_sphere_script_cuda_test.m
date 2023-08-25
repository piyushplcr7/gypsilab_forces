% Script for constant Magnetization
delete(gcp('nocreate'));
addpath(genpath("../../../"));
clear; clc; close all;
format long;

mu = 1;
vals = 5:13;
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
    bndmesh = getMeshSphere(N);
    
    % Mesh size
    hvals(i) = sqrt(mean(bndmesh.ndv,1));
    
    Gamma = dom(bndmesh,3);
    normals = Gamma.qudNrm;

    % Function to determine the interior of the magnet
    %interior = @(X) (vecnorm(X-T,2,2) < 1);
    
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
    [TnAJ,TdAJ,TnAM,TdAM] = solveMagnetProblemSimplified(Gamma,J,omega_src,mu,mu,M);

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

    %% Computing forces and torques from Equivalent current model

    % Force computation
    NED = fem(bndmesh,'NED'); 
    P1 = fem(bndmesh,'P1');
    P0 = fem(bndmesh,'P0');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0 = nxgrad(P1); 
    RWG = fem(bndmesh,'RWG');
    
    % Testing the condition number
%     testmat = integral(Gamma,DIV0,DIV0);
%     rcond(full(testmat))

    [X,W] = Gamma.qud;

    Mvals = M(X);
    Mxn = cross(Mvals,normals,2);

    Mdotn = dot(Mvals,normals,2);
    Mdotn_coeffs = proj(Mdotn,Gamma,P0);

    % Obtaining average B and H values on the magnet boundary

    % {B} = B.n + {Btan}
    Bn = reconstruct(TdAM+TdAJ,Gamma,curl(NED)).*normals;
    % nx curlAxn from the outside
    Btano = cross(normals,reconstruct(TnAJ+TnAM,Gamma,DIV0),2);
    Btani = Btano + Mxn;
    avgB = Bn + 0.5*(Btano + Btani);
    
    % Computing the integral of (Mxn)x{B}
    forces_mst(i,:) = sum(W.* cross(Mxn,avgB,2),1)

    Xcg = [4 0 0];
    r = X-Xcg;
    torques_mst(i,:) = sum(W.* cross(r,cross(Mxn,avgB,2),2),1)

    % Computing using BIE based formula
    [Vel1,DVel1] = getTransVelDVel([1 0 0]);
    [Vel2,DVel2] = getTransVelDVel([0 1 0]);
    [Vel3,DVel3] = getTransVelDVel([0 0 1]);

    fbem1 = PermanentMagnetShapeDerivativeBIEVPSimplified(Gamma,TnAJ+TnAM,TdAJ+TdAM,J,omega_src,Vel1);
    fbem2 = PermanentMagnetShapeDerivativeBIEVPSimplified(Gamma,TnAJ+TnAM,TdAJ+TdAM,J,omega_src,Vel2);
    fbem3 = PermanentMagnetShapeDerivativeBIEVPSimplified(Gamma,TnAJ+TnAM,TdAJ+TdAM,J,omega_src,Vel3);

    forces_bem(i,:) = [fbem1 fbem2 fbem3]

%     [Vel,DVel] = getPolyVelDVel(1,1,1,1);
    [Vel,DVel] = getCosVelDVel(1,1,1,1);

    tbem1 = PermanentMagnetShapeDerivativeBIEVP(Gamma,TnAJ_RWG+TnAM_RWG,TdAJ+TdAM,J,omega_src,mu,mu,M,Vel,DVel);
    
end
