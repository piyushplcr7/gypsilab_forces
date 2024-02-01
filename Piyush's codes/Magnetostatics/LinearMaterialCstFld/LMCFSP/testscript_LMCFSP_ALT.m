% Testscript LMCF SP

% New transmission problem script
delete(gcp('nocreate'));
addpath(genpath("../../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 4;
mu0 = 2;
vals = 5:12;
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
    
    bndmesh_i = getMeshSphere(N);

    % Mesh size
    hvals(i) = sqrt(mean(bndmesh_i.ndv,1));

    % Visualize with outer mesh translucent
    %customPlot(bndmesh_i,bndmesh_e);
    
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);
    normals_i = Gamma_i.qudNrm;
    normals_e = Gamma_e.qudNrm;
    
    %% Solving the transmission problem

    [psi_i,g_i,psi_e] = solveTPLMCFSP_ALT(bndmesh_i,bndmesh_e,mu,mu0,H0);

    % Inserting custom solution which holds for mu = mu0
    P1_i = fem(bndmesh_i,'P1');
    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');

    gradudotn = -reconstruct(psi_i,Gamma_i,P0_i);
    sgradu = reconstruct(g_i,Gamma_i,grad(P1_i));

    Ht = vecnorm(sgradu,2,2);
    Bn = mu0 * gradudotn;



    %% Checking solution for mu = mu0 using the field
    a = 0; b = 1; c = 1; alpha = 0; kappa = 3;
    idx = a + kappa * b + kappa^2 * c + kappa^3 * alpha + 1
    [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);

    testsdmst(i) = ShapeDervTpVol(Gamma_i,Bn,Ht,mu0,mu,Vel)

    testsdbem(i) = SdBEMLMCFSP_ALT(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel,DVel,mu0,mu,H0)
    

end
