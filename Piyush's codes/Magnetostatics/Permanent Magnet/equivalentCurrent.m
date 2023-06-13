% Script for constant Magnetization

addpath(genpath("../../../"));
clear; clc; close all;
format long;

mu = 1;
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
    interior = @(X) (vecnorm(X-T,2,2) < 1);
    
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
    
    [X,W] = Gamma.qud;

    Mvals = M(X);
    Mxn = cross(Mvals,normals,2);
    Mxn_coeffs = proj(Mxn,Gamma,DIV0);

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
    force = sum(W.* cross(Mxn,avgB,2),1)

    Xcg = [4 0 0];
    r = X-Xcg;
    torque = sum(W.* cross(r,cross(Mxn,avgB,2),2),1)

end
