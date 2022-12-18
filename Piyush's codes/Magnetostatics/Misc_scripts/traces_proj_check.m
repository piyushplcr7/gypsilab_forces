% Script checking convergence of the projected traces obtained from a
% source current

addpath(genpath("../../"));

clear; clc;

for i = 3:11
    N = 2^i
    
    % Cube size and position
    L = 5*[1 1 1];
    T = [10 0 0];
    
    % Solution domain
    %bndmesh = bndmeshCubeTranslated(N,L,T);
    bndmesh = mshSphere(N,2);
    bndmesh = bndmesh.translate(T);
    Gamma = dom(bndmesh,3);

    % BEM spaces
    NED = fem(bndmesh,'NED'); % curl conforming -> Dirichlet trace
    P1 = fem(bndmesh,'P1');
    DIV0 = nxgrad(P1); % div conforming, div0 -> Neumann trace
    DIV = fem(bndmesh,'RWG');

    %% Generating synthetic traces from a source
    N_src = N;
    R0 = 2;
    r0 = .5;
    [J,mesh_src] = get_torus_source(N_src,R0,r0);
    omega_src = dom(mesh_src,3);

    % Evaluation points on Gamma
    [X,W,elt2qud] = Gamma.qud;

    % Tangential vector fields
    A = compute_vecpot(J,omega_src,X);
    curlA = compute_vecpot_curl(J,omega_src,X);

    normals = Gamma.qudNrm;

    TDA = A - dot(A,normals,2).*normals;
    TNA = cross(normals,curlA,2);

    %% Projecting traces and comparing
    TDA_NED_coeffs = proj(TDA,Gamma,NED);
    TDA_DIV_coeffs = proj(TDA,Gamma,DIV);
    %TDA_DIV0_coeffs = proj(TDA,Gamma,DIV0);

    TDA_NED_recon = reconstruct(TDA_NED_coeffs,Gamma,NED);
    TDA_DIV_recon = reconstruct(TDA_DIV_coeffs,Gamma,DIV);
    %TDA_DIV0_recon = reconstruct(TDA_DIV0_coeffs,Gamma,DIV0);

    %quiver3wrapper(X,TDA_NED_recon,'blue');
    %hold on;
    %quiver3wrapper(X,TDA,'red');
    %quiver3wrapper(X,TDA_DIV_recon,'green');

    L2err_NED = sum(W.*(vecnorm(TDA-TDA_NED_recon,2,2).^2),1)
    L2err_DIV = sum(W.*(vecnorm(TDA-TDA_DIV_recon,2,2).^2),1)
    
    % Unknown error metric lol
    uerrned = norm(TDA-TDA_NED_recon)
    uerrdiv = norm(TDA-TDA_DIV_recon)
    %norm(TDA-evalTDA)%/norm(TDA)
    %norm(TNA-evalTNA)%/norm(TNA)
end