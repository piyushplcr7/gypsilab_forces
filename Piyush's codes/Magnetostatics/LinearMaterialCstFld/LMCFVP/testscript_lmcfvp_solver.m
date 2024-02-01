% New transmission problem script
delete(gcp('nocreate'));
addpath(genpath("../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 1;
mu0 = mu;
vals = 5:9;
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
forces_bem = forces_mst;
torques_mst = forces_mst;
torques_bem = forces_mst;
hvals = 0*vals;

for i = 1:Nvals
    N = 2^vals(i);
    disp(N);
    %% SOLUTION DOMAIN
    % Cube size and position
    L = 2*[1 1 1];
    T = [1 0.5 2];
    
    % Bounding box
    bndmesh_e = mshSphere(N,5);
    
    bndmesh_i = mshSphere(N,1);


%     bndmesh_i = mshSphere(N,1);
%     bndmesh_i = bndmesh_i.translate([2 0 0]);

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
    
    %% Solving the transmission problem
    
    B_0 = [1 0 0];
    A_exact = @(X) 0.5 * cross(repmat(B_0,size(X,1),1),X,2);

    % A_exact = @(X) ones(size(X,1),1) * B_0;

    NED_e = fem(bndmesh_e,'NED');
    [X_e,~] = Gamma_e.qud;
    normals_e = Gamma_e.qudNrm;
    g_e_vals = cross(normals_e, cross(A_exact(X_e),normals_e,2) ,2);
    g_e = proj(g_e_vals,Gamma_e,NED_e);
    
    % These are traces from the exterior
    [Psi_i,g_i,Psi_e] = solveTPLMCFVP_ALT(bndmesh_i,bndmesh_e,mu,mu0,g_e);

    %% Exact computations
    NED_i = fem(bndmesh_i,'NED'); 
    P1_i = fem(bndmesh_i,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0_i = nxgrad(P1_i); 
    
    % Bn
    curlAdotn = reconstruct(g_i,Gamma_i,NED_i.curl);
    curlAtan = -cross(normals_i,reconstruct(Psi_i,Gamma_i,DIV0_i),2);

    curlA_reconstructed = curlAdotn.*normals_i + curlAtan;

    % error_curlA = curlA_reconstructed - B_0;
    error_curlA = curlA_reconstructed;

    errnorm = vecnorm(error_curlA,2,2);

    plot(errnorm);

    B0_extended = repmat(B_0,size(normals_i,1),1);
    [X_i,~] = Gamma_i.qud;

    figure;
    quiver3wrapper(X_i,B0_extended,'red');
    axis equal;
    hold on;
    quiver3wrapper(X_i,curlA_reconstructed,'blue');
    
    disp('abc');

end
