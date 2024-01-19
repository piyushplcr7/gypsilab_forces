delete(gcp('nocreate'));
addpath(genpath("../../../"));
format long;
mu0 = 1;
vals = 10:12;
Nvals = size(vals,2);

forces_mst = zeros(Nvals,3);
forces_bem = forces_mst;
torques_mst = forces_mst;
torques_bem = forces_mst;
hvals = 0*vals;
testsdmst = hvals;
testsdbem = hvals;

for i = 1:Nvals
    N = 2^vals(i);
    disp(N);
    %% SOLUTION DOMAIN
    r_i = 2;
    r_o = 4;

    bndmesh_i = mshSphere(N,r_i);

    % Bounding box
    bndmesh_e = mshSphere(N,r_o);
    
    Gamma_i = dom(bndmesh_i,3);
    normals_i = Gamma_i.qudNrm;
    
    %% Solving the transmission problem
    jump = -1/r_i^2;
    P0_i = fem(bndmesh_i,'P0');
    P1_i = fem(bndmesh_i,'P1');
    [X_i,W_i] = Gamma_i.qud;

    jump_extended = repelem(jump,size(X_i,1),1);
    jump_coeffs = proj(jump_extended,Gamma_i,P0_i);
    
    [psi_i,g_i,psi_e] = solveTpPMCFSP_Jump(bndmesh_i,bndmesh_e,jump_coeffs);

    %% Comparing with exact solution

    g_exact_i = 1/r_i - 1/r_o;
    psi_exact_i = 1/r_i^2;

    g_i_reconstructed = reconstruct(g_i,Gamma_i,P1_i);
    psi_i_recostructed = reconstruct(psi_i,Gamma_i,P0_i);
    
   
end