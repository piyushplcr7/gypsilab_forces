% Testscript LMCF SP

% New transmission problem script
delete(gcp('nocreate'));
addpath(genpath("../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 4;
mu0 = 2;
vals = 5:9;
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
forces_mst_recon = forces_mst;
forces_bem = forces_mst;
torques_mst = forces_mst;
torques_mst_recon = torques_mst;
torques_bem = forces_mst;
hvals = 0*vals;
testsdbem = hvals;
testsdmst = testsdbem;

%rng(32);
H0 = [10 3 1];

for i = 1:Nvals
    
    N = 2^vals(i);
    disp(N);
    %% SOLUTION DOMAIN
    % Cube size and position
    L = 2*[1 1 1];
    T = [1 0.5 2];
    
    % Bounding box
    bndmesh_e = mshSphere(N,9);
    bndmesh_e = bndmesh_e.translate([2 2 2]);
    
    mesh = mshCube(N,L);
    mesh = mesh.translate(T);
    %mesh = mesh.sub(1);
    bndmesh_i = mesh.bnd;

    bndmesh_i = getMeshSphere(N);

    % bndmesh_i = mshSphere(N,1);
    % bndmesh_i = bndmesh_i.translate([2 0 0]);
%     cutoff_rad = 3.95;

%     assert(max(vecnorm(bndmesh_i.vtx,2,2))<cutoff_rad);
%     assert(min(vecnorm(bndmesh_e.vtx,2,2))>cutoff_rad);
    % Mesh size
    hvals(i) = sqrt(mean(bndmesh_i.ndv,1));

    % Visualize with outer mesh translucent
    %customPlot(bndmesh_i,bndmesh_e);
    
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);
    normals_i = Gamma_i.qudNrm;
    
    %% Solving the transmission problem

%     H0 = [1 0 0];%[1 1 1]/sqrt(3);

    % Match for H0 = 0 1 0
    % [psi_i,g_i,psi_e,psi_I_recon,g_i_recon] = solveTPLMCFSP(bndmesh_i,bndmesh_e,mu,mu0,H0);
    [psi_i,g_i,psi_e] = solveTPLMCFSP(bndmesh_i,bndmesh_e,mu,mu0,H0);
    %% Computing the force using MST formula
    % Reconstructing Bn and Ht
%     P1_i = fem(bndmesh_i,'P1');
%     P0_i = fem(bndmesh_i,'P0');
%     H0extended = repmat(H0,size(normals_i,1),1);
%     H0t = cross(normals_i,cross(H0extended,normals_i,2),2);
%     Ht = reconstruct(g_i,Gamma_i,P1_i.grad) + H0t;
% 
%     Bn = -mu0 * reconstruct(psi_i,Gamma_i,P0_i) + mu0 * dot(H0extended,normals_i,2);
% 
%     forces_mst(i,:) = ForceMstTP(Gamma_i,Bn,vecnorm(Ht,2,2),mu0,mu)
% 
%     Xcg = [4 0 0];
%     torques_mst(i,:) = TorqueMstTP(Gamma_i,Bn,vecnorm(Ht,2,2),mu0,mu,Xcg)

    %% Computing the force using MST formula (reconstructed traces)
    % Reconstructing Bn and Ht
    % Ht = reconstruct(g_i_recon,Gamma_i,P1_i.grad) + H0t;
    % 
    % Bn = mu * reconstruct(psi_I_recon,Gamma_i,P0_i) + mu0 * dot(H0extended,normals_i,2);
    % 
    % forces_mst_recon(i,:) = ForceMstTP(Gamma_i,Bn,vecnorm(Ht,2,2),mu0,mu)
    % 
    % Xcg = [4 0 0];
    % torques_mst_recon(i,:) = TorqueMstTP(Gamma_i,Bn,vecnorm(Ht,2,2),mu0,mu,Xcg)

    %% Computing force using BEM formula

%     [Vel1,DVel1] = getTransVelDVelCutoff([1 0 0],cutoff_rad);
%     [Vel2,DVel2] = getTransVelDVelCutoff([0 1 0],cutoff_rad);
%     [Vel3,DVel3] = getTransVelDVelCutoff([0 0 1],cutoff_rad);
% 
%     f1 = SdBEMLMCFSP_ConstVEL(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel1,DVel1,mu0,mu,H0);
%     f2 = SdBEMLMCFSP_ConstVEL(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel2,DVel2,mu0,mu,H0);
%     f3 = SdBEMLMCFSP_ConstVEL(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel3,DVel3,mu0,mu,H0);

%     f1 = SdBEMLMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel1,DVel1,mu0,mu,H0);
%     f2 = SdBEMLMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel2,DVel2,mu0,mu,H0);
%     f3 = SdBEMLMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel3,DVel3,mu0,mu,H0);

%     forces_bem(i,:) = [f1 f2 f3]

    % Reconstructing Bn and Ht
        P1_i = fem(bndmesh_i,'P1');
        P0_i = fem(bndmesh_i,'P0');
        H0extended = repmat(H0,size(normals_i,1),1);
        H0t = cross(normals_i,cross(H0extended,normals_i,2),2);
        Ht = reconstruct(g_i,Gamma_i,P1_i.grad) + H0t;
        Ht = vecnorm(Ht,2,2);
    
        Bn = -mu0 * reconstruct(psi_i,Gamma_i,P0_i) + mu0 * dot(H0extended,normals_i,2);
        jump_mu_inv = 1/mu0 - 1/mu;
        jump_mu = mu0 - mu;
        [X_i,W_i] = Gamma_i.qud;
        fdensity = 0.5 * ((Bn).^2*jump_mu_inv - (Ht).^2*jump_mu).* normals_i;

    a = 1; b = 2; c = 2; alpha = 1; kappa = 3;
    idx = a + kappa * b + kappa^2 * c + kappa^3 * alpha + 1;
    [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);
    Vels = Vel(X_i);

    testsdmst(i) = sum(W_i.*dot(fdensity,Vels,2),1)
    testsdbem(i) = SdBEMLMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel,DVel,mu0,mu,H0)

%     [Velr1,DVelr1] = getRotVelDVelCutoff([1 0 0],Xcg,cutoff_rad);
%     [Velr2,DVelr2] = getRotVelDVelCutoff([0 1 0],Xcg,cutoff_rad);
%     [Velr3,DVelr3] = getRotVelDVelCutoff([0 0 1],Xcg,cutoff_rad);
% 
%     t1 = SdBEMLMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Velr1,DVelr1,mu0,mu,H0);
%     t2 = SdBEMLMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Velr2,DVelr2,mu0,mu,H0);
%     t3 = SdBEMLMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Velr3,DVelr3,mu0,mu,H0);
% 
%     torques_bem(i,:) = [t1 t2 t3]


end
