% Solving the superconductor scalar potential Dirichlet Problem 

function [Tdu,Tnu] = solve_superconductor_scalar_potential(bndmesh,J,mesh_src)

    % BEM Spaces
    P1 = fem(bndmesh,'P1');
    P0 = fem(bndmesh,'P0');
    % Surface curl of P1
    SCurlP1 = nxgrad(P1);

    % DOM objects
    Gamma = dom(bndmesh,3);
    omega_src = dom(mesh_src,3);
    [X,W] = Gamma.qud;
    normals = Gamma.qudNrm;

    % Assembling the LHS Galerkin matrix
    Wmat = single_layer(Gamma,SCurlP1,SCurlP1);

    % Assembling the RHS
    % Computing HJ at the quadrature points on Gamma
    HJ = compute_vecpot_curl(J,omega_src,X);
    % Taking normal trace on Gamma
    HJn = dot(HJ,normals,2);
    % Projecting to P0
    HJn_coeffs = proj(HJn,Gamma,P0);

    M = mass_matrix(Gamma,P1,P0);
    Kmat = double_layer_laplace(Gamma,P0,P1);

    rhs = (0.5*M+Kmat')*HJn_coeffs;

    % vector to enforce zero mean constraint
    vec = integral(Gamma,P1);
    blockmat = [Wmat vec;
                vec' 0];
    rhsmodif = [rhs;0];

    solextended = blockmat\rhsmodif;

    Tdu = solextended(1:end-1); 
    Tnu = -HJn_coeffs;

end