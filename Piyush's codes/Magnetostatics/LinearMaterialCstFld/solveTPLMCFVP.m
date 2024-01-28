% Solver linear material constant field vector potential
% Solve Transmission Problem. Returns the exterior traces 

function [Psi_i,g_i,Psi_e] = solveTPLMCFVP(bndmesh_i,bndmesh_e,mu,mu_0,B_0)
    %% BEM Spaces
    Gamma_i = dom(bndmesh_i,3);
    [X_i,~] = Gamma_i.qud;

    Gamma_e = dom(bndmesh_e,3);
    [X_e,~] = Gamma_e.qud;

    % BEM spaces
    % curl conforming -> Dirichlet trace
    NED_i = fem(bndmesh_i,'NED'); 
    P1_i = fem(bndmesh_i,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0_i = nxgrad(P1_i); 
    % Kernel of the surface curl operator
    Ker_curl_i = grad(P1_i); 
    % Div conforming space 
    DIV_i = fem(bndmesh_i,'RWG');

    NED_e = fem(bndmesh_e,'NED');
    P1_e = fem(bndmesh_e,'P1');
    DIV0_e = nxgrad(P1_e);
    DIV_e = fem(bndmesh_e,'RWG');

    %% Galerkin matrices 
    % For operator A
    Aii = single_layer(Gamma_i,DIV0_i,DIV0_i);
    % For operator C
    Cii = double_layer_magnetostatics(Gamma_i,DIV0_i,DIV_i);

    %% Linear System

    blockopr = LMCFVP_systemMat(bndmesh_i,bndmesh_e,mu,mu_0);

    %% Computing the RHS

    % Computing B_0 X n at quadrature points on Gamma_i
    normals_i = Gamma_i.qudNrm;

    B_0xn = cross(repmat([B_0(1) B_0(2) B_0(3)],size(normals_i,1),1),normals_i,2);
    B_0xn_coeffs = proj(B_0xn,Gamma_i,DIV0_i);
    % Jump
    jumpMuInv = 1/mu_0 - 1/mu;

    % rhs1
    rhs1 = mu * jumpMuInv * Aii * B_0xn_coeffs;

    % rhs2
    massmat = mass_matrix(Gamma_i,NED_i,DIV0_i);
    rhs2 = -mu_0 * jumpMuInv * Cii' * B_0xn_coeffs + mu_0 * jumpMuInv/2 * massmat * B_0xn_coeffs;

    rhs = [rhs1; rhs2; zeros(P1_e.ndof,1) ;0; zeros(P1_i.ndof,1); 0; 0];
    sol = blockopr\rhs;

    Psi_i = sol(1:P1_i.ndof);
    g_i = sol(P1_i.ndof+1:P1_i.ndof+DIV_i.ndof);
    Psi_e = sol(P1_i.ndof+DIV_i.ndof+1: P1_i.ndof+DIV_i.ndof + P1_e.ndof);

    %% Verification of the solution via BIEs

end