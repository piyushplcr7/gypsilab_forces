% Testscript LMCF SP

% New transmission problem script
delete(gcp('nocreate'));
addpath(genpath("../../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 24;
mu0 = 2;
vals = 8:9;
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
    
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);
    normals_i = Gamma_i.qudNrm;
    normals_e = Gamma_e.qudNrm;

    P1_i = fem(bndmesh_i,'P1');
    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');
    
    %% LMCFSP alt solution and reconstruction

    [psi_i_ALT,g_i_ALT,psi_e_ALT] = solveTPLMCFSP_ALT(bndmesh_i,bndmesh_e,mu,mu0,H0);

    gradudotn_ALT = -reconstruct(psi_i_ALT,Gamma_i,P0_i);
    sgradu_ALT = reconstruct(g_i_ALT,Gamma_i,grad(P1_i));

    Ht = vecnorm(sgradu_ALT,2,2);
    Bn = mu0 * gradudotn_ALT;

    Halt = Bn/mu0.*normals_i + sgradu_ALT;

    %% LMCFSP solution and reconstruction


    [psi_i,g_i,psi_e] = solveTPLMCFSP(bndmesh_i,bndmesh_e,mu,mu0,H0); 

    H0extended = repmat(H0,size(normals_i,1),1);
    H0t = cross(normals_i,cross(H0extended,normals_i,2),2);
    Htot_t = reconstruct(g_i,Gamma_i,P1_i.grad) + H0t;
    % Htot_t = vecnorm(Htot_t,2,2);

    Bntot = -mu0 * reconstruct(psi_i,Gamma_i,P0_i) + mu0 * dot(H0extended,normals_i,2);

    Hreg = Htot_t + Bntot/mu0.*normals_i; 
    
    %% Visualizing
    [X_i,W_i] = Gamma_i.qud;
    quiver3wrapper(X_i,Halt,'red');
    hold on;
    quiver3wrapper(X_i,Hreg,'blue');


    %% Comparing the solutions
    

    l2err_bn = sum(W_i.*(Bn-Bntot).^2,1)
    l2err_ht = sum(W_i.*vecnorm(sgradu_ALT-Htot_t,2,2).^2,1)

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
