delete(gcp('nocreate'));
addpath(genpath("../../../"));
format long;
mu0 = 1;
vals = 5:12;
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
    bndmesh_i = getMeshSphere(N);

    % Bounding box
    bndmesh_e = mshSphere(N,9);
    bndmesh_e = bndmesh_e.translate([2 2 2]);
    assert(max(vecnorm(bndmesh_i.vtx-[2 2 2],2,2))<9);
    % Mesh size
    hvals(i) = sqrt(mean(bndmesh_i.ndv,1));
    
    Gamma_i = dom(bndmesh_i,3);
    normals_i = Gamma_i.qudNrm;
    
    %% Solving the transmission problem
    M = @(X) ones(size(X,1),1)  * [10 0 0];
    H0 = [1 0 0];
    [psi_i,g_i,psi_e] = solveTpPMCFSP(bndmesh_i,bndmesh_e,mu0,M);

    %% Computing the force using equivalent charge model
    P1_i = fem(bndmesh_i,'P1');
    P0_i = fem(bndmesh_i,'P0');
    H0extended = repmat(H0,size(normals_i,1),1);

    [X_i,W_i] = Gamma_i.qud;

    Mvals = M(X_i);
    Mdotn = dot(Mvals,normals_i,2);

    psi_i_vals = reconstruct(psi_i,Gamma_i,P0_i);
    gradg_i_vals = reconstruct(g_i,Gamma_i,P1_i.grad);

    avgHtot = (-psi_i_vals-0.5 * Mdotn).*normals_i + gradg_i_vals + H0extended;

    fdensity = avgHtot .* Mdotn;

    %%

    a = 2;
    b = 2;
    c = 2;
    alpha = 0;
    kappa = 3;
    idx = a + kappa * b + kappa^2 * c + kappa^3 * alpha + 1;

    [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);
    Vels = Vel(X_i);

    testsdmst(i) = mu0 * sum(W_i.*dot(fdensity,Vels,2),1)

    testsdbem(i) = SdBEMPMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel,DVel,mu0,H0,M)
    
   
end