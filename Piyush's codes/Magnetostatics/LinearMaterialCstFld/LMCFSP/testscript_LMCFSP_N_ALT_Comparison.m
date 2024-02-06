% Testscript LMCF SP

% New transmission problem script
delete(gcp('nocreate'));
addpath(genpath("../../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 4;
mu0 = 2;
vals = 5:10;
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
forces_mst_recon = forces_mst;
forces_bem = forces_mst;
torques_mst = forces_mst;
torques_mst_recon = torques_mst;
torques_bem = forces_mst;
hvals = 0*vals;

H0 = [10 3 1];

for i = 1:Nvals
    N = 2^vals(i);
    disp(N);
    %% SOLUTION DOMAIN
    
    % Bounding box
    bndmesh_e = mshSphere(N,9);
    bndmesh_e = bndmesh_e.translate([2 2 2]);
    
    bndmesh_i = getMeshCube(N);

    % Mesh size
    hvals(i) = sqrt(mean(bndmesh_i.ndv,1));

    % Visualize with outer mesh translucent
    %customPlot(bndmesh_i,bndmesh_e);
    order = 7;
    Gamma_i = dom(bndmesh_i,order);
    Gamma_e = dom(bndmesh_e,order);
    normals_i = Gamma_i.qudNrm;
    normals_e = Gamma_e.qudNrm;

    P1_i = fem(bndmesh_i,'P1');
    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');
    
    %% LMCFSP alt solution and reconstruction

    [psi_i_ALT,g_i_ALT,psi_e_ALT] = solveTPLMCFSP_ALT(bndmesh_i,bndmesh_e,mu,mu0,H0);

    gradudotn_ALT = -reconstruct(psi_i_ALT,Gamma_i,P0_i);
    sgradu_ALT = reconstruct(g_i_ALT,Gamma_i,grad(P1_i));

    Ht_ALT = vecnorm(sgradu_ALT,2,2);
    Bn_ALT = mu0 * gradudotn_ALT;

    Halt = Bn_ALT/mu0.*normals_i + sgradu_ALT;

    %% LMCFSP solution and reconstruction


    [psi_i,g_i,psi_e] = solveTPLMCFSP(bndmesh_i,bndmesh_e,mu,mu0,H0); 

    H0extended = repmat(H0,size(normals_i,1),1);
    H0t = cross(normals_i,cross(H0extended,normals_i,2),2);
    Htot_t = reconstruct(g_i,Gamma_i,P1_i.grad) + H0t;

    Bntot = -mu0 * reconstruct(psi_i,Gamma_i,P0_i) + mu0 * dot(H0extended,normals_i,2);

    H = Htot_t + Bntot/mu0.*normals_i; 
    
    
    %% Visualizing
    [X_i,W_i] = Gamma_i.qud;
    % quiver3wrapper(X_i,Halt,'red');
    % hold on;
    % quiver3wrapper(X_i,H,'blue');


    %% Comparing the solutions (L2 norm of difference)
    

    % l2err_bn = sum(W_i.*(Bn_ALT-Bntot).^2,1)
    % l2err_ht = sum(W_i.*vecnorm(sgradu_ALT-Htot_t,2,2).^2,1)

    %% Comparing the MST shape derivatives

    jump_mu_inv = 1/mu0 - 1/mu;
    jump_mu = mu0 - mu;

    fdensity_ALT = 0.5 * ((Bn_ALT).^2*jump_mu_inv - (Ht_ALT).^2*jump_mu).* normals_i;
    fdensity = 0.5 * ((Bntot).^2*jump_mu_inv - (vecnorm(Htot_t,2,2)).^2*jump_mu).* normals_i;

    [Vel,DVel] = getTransVelDVel([1 0 0]);
    Vels = Vel(X_i);

    testsd_mst(i) = sum(W_i.*dot(fdensity,Vels,2),1)
    testsd_mst_ALT(i) = sum(W_i.*dot(fdensity_ALT,Vels,2),1)
    

    %% Checking solution for mu = mu0 using the field
    % a = 0; b = 1; c = 1; alpha = 0; kappa = 3;
    % idx = a + kappa * b + kappa^2 * c + kappa^3 * alpha + 1
    % [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);
    % 
    % [Vel,DVel] = getTransVelDVel([1 0 0]);
    % 
    % testsdmst(i) = ShapeDervTpVol(Gamma_i,Bn,Ht,mu0,mu,Vel)

    % testsdbem_constvel(i) = SdBEMLMCFSP_ConstVEL_ALT(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel,DVel,mu0,mu,H0)

    % testsdbem(i) = SdBEMLMCFSP_ALT(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel,DVel,mu0,mu,H0)
    

end
