% Testscript LMCF SP

% New transmission problem script
delete(gcp('nocreate'));
addpath(genpath("../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 1;
mu0 = 1;
vals = 5:9;
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
forces_mst_recon = forces_mst;
forces_bem = forces_mst;
torques_mst = forces_mst;
torques_mst_recon = torques_mst;
torques_bem = forces_mst;
hvals = 0*vals;

rng(32);
H0 = rand(1,3);

for i = 1:Nvals
    N = 2^vals(i);
    disp(N);
    %% SOLUTION DOMAIN
    % Cube size and position
    L = 2*[1 1 1];
    T = [1 0.5 2];
    % T = [1 1.5 0.5];
    
    % Bounding box
    bndmesh_e = mshSphere(2*N,5);
    
    mesh = mshCube(N,L);
    mesh = mesh.translate(T);
    bndmesh_i = mesh.bnd;

    % bndmesh_i = mshSphere(N,1);
    % bndmesh_i = bndmesh_i.translate([2 0 0]);
    cutoff_rad = 3.95;

    assert(max(vecnorm(bndmesh_i.vtx,2,2))<cutoff_rad);
    assert(min(vecnorm(bndmesh_e.vtx,2,2))>cutoff_rad);
    % Mesh size
    hvals(i) = sqrt(mean(bndmesh_i.ndv,1));

    % Visualize with outer mesh translucent
    %customPlot(bndmesh_i,bndmesh_e);
    
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);
    normals_i = Gamma_i.qudNrm;
    normals_e = Gamma_e.qudNrm;
    
    %% Solving the transmission problem
    %H0 = [1 1 1]/sqrt(3);

    % Match for H0 = 0 1 0
    [psi_i,g_i,psi_e] = solveTPLMCFSP_ALT(bndmesh_i,bndmesh_e,mu,mu0,H0);

    % Inserting custom solution which holds for mu = mu0
    P1_i = fem(bndmesh_i,'P1');
    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');
     expl_psi_i_vals = -normals_i * H0';
     psi_i = proj(expl_psi_i_vals,Gamma_i,P0_i);
     [X_i,~] = Gamma_i.qud;
     expl_g_i_vals = X_i * H0';
     g_i = proj(expl_g_i_vals,Gamma_i,P1_i);
     expl_psi_e_vals = normals_e * H0';
     psi_e = proj(expl_psi_e_vals,Gamma_e,P0_e);

    %% Checking something 
    % [X_i,W_i] = Gamma_i.qud;
    % psi_i_vals = reconstruct(psi_i,Gamma_i,P0_i);
    % int_psi_i = sum(W_i.*psi_i_vals,1);
    % int_Gamma_i = sum(W_i,1);
    % altint_Gamma_i = sum(bndmesh_i.ndv,1);
    % disp(int_psi_i/altint_Gamma_i);
    
    %% Computing the energy
    % normals_e = Gamma_e.qudNrm;
    % H0dotn_e = normals_e * H0';
    % P0_e = fem(bndmesh_e,'P0');
    % P1_e = fem(bndmesh_e,'P1');
    % [X_e,~] = Gamma_e.qud;
    % g_e_vals = X_e * H0';
    % g_e_coeffs = proj(g_e_vals,Gamma_e,P1_e);
    % energy = mu0/2 * psi_e' * mass_matrix(Gamma_e,P0_e,P1_e) * g_e_coeffs;
    % H0dotn_e_coeffs = proj(H0dotn_e,Gamma_e,P0_e);
    % energy_exact = mu0/2 * H0dotn_e_coeffs' * mass_matrix(Gamma_e,P0_e,P1_e) * g_e_coeffs;
    % disp(energy);
    % disp(energy_exact);
    %% Computing the force using MST formula

    Ht = reconstruct(g_i,Gamma_i,P1_i.grad);

    Bn = -mu0 * reconstruct(psi_i,Gamma_i,P0_i);

    forces_mst(i,:) = ForceMstTP(Gamma_i,Bn,vecnorm(Ht,2,2),mu0,mu)
    % 
    Xcg = [4 0 0];
    torques_mst(i,:) = TorqueMstTP(Gamma_i,Bn,vecnorm(Ht,2,2),mu0,mu,Xcg)

    %% Computing force using BEM formula

    [Vel1,DVel1] = getTransVelDVelCutoff([1 0 0],cutoff_rad);
    [Vel2,DVel2] = getTransVelDVelCutoff([0 1 0],cutoff_rad);
    [Vel3,DVel3] = getTransVelDVelCutoff([0 0 1],cutoff_rad);

    f1 = SdBEMLMCFSP_ConstVEL_ALT(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel1,DVel1,mu0,mu,H0);
    f2 = SdBEMLMCFSP_ConstVEL_ALT(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel2,DVel2,mu0,mu,H0);
    f3 = SdBEMLMCFSP_ConstVEL_ALT(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel3,DVel3,mu0,mu,H0);

    forces_bem(i,:) = [f1 f2 f3]

    [Velr1,DVelr1] = getRotVelDVelCutoff([1 0 0],Xcg,cutoff_rad);
    [Velr2,DVelr2] = getRotVelDVelCutoff([0 1 0],Xcg,cutoff_rad);
    [Velr3,DVelr3] = getRotVelDVelCutoff([0 0 1],Xcg,cutoff_rad);

    t1 = SdBEMLMCFSP_ALT(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Velr1,DVelr1,mu0,mu,H0);
    t2 = SdBEMLMCFSP_ALT(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Velr2,DVelr2,mu0,mu,H0);
    t3 = SdBEMLMCFSP_ALT(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Velr3,DVelr3,mu0,mu,H0);

    torques_bem(i,:) = [t1 t2 t3]
    

end
