% Testscript LMCF SP

% New transmission problem script
delete(gcp('nocreate'));
addpath(genpath("../../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 2;
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

    bndmesh_i = getMeshSphere(N);

    % Mesh size
    hvals(i) = sqrt(mean(bndmesh_i.ndv,1));

    % Visualize with outer mesh translucent
    %customPlot(bndmesh_i,bndmesh_e);
    
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);
    normals_i = Gamma_i.qudNrm;
    
    %% Solving the transmission problem

    [psi_i,g_i,psi_e] = solveTPLMCFSP_ALT(bndmesh_i,bndmesh_e,mu,mu0,H0); 
    

    %% Computing the force using MST formula (reconstructed traces

    % Reconstructing Bn and Ht
    P1_i = fem(bndmesh_i,'P1');
    P0_i = fem(bndmesh_i,'P0');
    
    Ht = reconstruct(g_i,Gamma_i,P1_i.grad);

    Hn = -reconstruct(psi_i,Gamma_i,P0_i);

    Htot = Hn.*normals_i + Ht;

    % plot(bndmesh_i);
    hold on;

    H0extended = repmat(H0,size(normals_i,1),1);

    [X_i,~] = Gamma_i.qud;
    quiver3wrapper(X_i,H0extended,'red');

    hold on;
    quiver3wrapper(X_i,Htot,'blue');


    disp('yo')


end
