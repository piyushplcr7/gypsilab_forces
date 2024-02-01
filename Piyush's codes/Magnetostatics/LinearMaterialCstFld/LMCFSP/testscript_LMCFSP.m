% Testscript LMCF SP

% New transmission problem script
delete(gcp('nocreate'));
addpath(genpath("../../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 4;
mu0 = 2;
vals = 7:12;
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
% H0 = [1 0 0];

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
    
    %% Solving the transmission problem

    [psi_i,g_i,psi_e] = solveTPLMCFSP(bndmesh_i,bndmesh_e,mu,mu0,H0); 
    

    %% Computing the force using MST formula (reconstructed traces

    % Reconstructing Bn and Ht
    P1_i = fem(bndmesh_i,'P1');
    P0_i = fem(bndmesh_i,'P0');
    H0extended = repmat(H0,size(normals_i,1),1);
    H0t = cross(normals_i,cross(H0extended,normals_i,2),2);
    Htot_t = reconstruct(g_i,Gamma_i,P1_i.grad) + H0t;
    Htot_t = vecnorm(Htot_t,2,2);

    Bntot = -mu0 * reconstruct(psi_i,Gamma_i,P0_i) + mu0 * dot(H0extended,normals_i,2);
    jump_nu = 1/mu0 - 1/mu;
    jump_mu = mu0 - mu;
    [X_i,W_i] = Gamma_i.qud;
    fdensity = 0.5 * (jump_nu * (Bntot).^2 - jump_mu * (Htot_t).^2).* normals_i;

    a = 0; b = 1; c = 1; alpha = 0; kappa = 3;
    idx = a + kappa * b + kappa^2 * c + kappa^3 * alpha + 1
    [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);

    [Vel,DVel] = getRotVelDVel([1 0 0],[5 5 3]);
    % [Vel,DVel] = getRotVelDVel([1 0 0],[4 0 0]);
    % [Vel,DVel] = getTransVelDVel([1 0 0]);
    Vels = Vel(X_i);

    % testsdmst(i) = sum(W_i.*dot(fdensity,Vels,2),1)
    testsdmst(i) = ShapeDervTpVol(Gamma_i,Bntot,Htot_t,mu0,mu,Vel)
    testsdbem(i) = SdBEMLMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel,DVel,mu0,mu,H0)



end
