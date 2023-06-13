% Solving transmission problem for permanent magnet using scalar potential
% This script assumes that divM = 0

function [psi,g] = solveMagnetProblemSimplifiedSP(Gamma,mu,mu0,J,omega_src,M)
    % BEM spaces
    bndmesh = Gamma.msh;
    P1 = fem(bndmesh,'P1');
    P0 = fem(bndmesh,'P0');

    [X,~] = Gamma.qud;
    normals = Gamma.qudNrm;

    V = single_layer(Gamma,P0,P0);
    K = double_layer_laplace(Gamma,P0,P1);
    W = single_layer(Gamma,P1.nxgrad,P1.nxgrad);
    % Vector to enforce zero P1 mean constraint
    Vec = integral(Gamma,P1);

    jumpMu = mu0-mu;

    HJ = compute_vecpot_curl(J,omega_src,X);
    HJn = dot(normals,HJ,2);
    HJn_coeffs = proj(HJn,Gamma,P0);

    % Evaluating M at the boundary
    Mbnd = M(X);
    Mdotn = dot(normals,Mbnd,2);
    Mdotn_coeffs = proj(Mdotn,Gamma,P0);

    % alpha = mu M.n - [mu] Hj.n
    alpha_coeffs = mu * Mdotn_coeffs - jumpMu * HJn_coeffs;

    %% Solving linear system
    NP1 = P1.ndof;
    NP0 = P0.ndof;
    % LHS matrix
    blockmat = [(1+mu0/mu)*V -2*K zeros(NP0,1);
                2*K' (1+mu/mu0)*W Vec;
                zeros(1,NP0) Vec' 0];

    % RHS
    rhs1 = 1/mu * V * alpha_coeffs;
    M_P1_P0 = mass_matrix(Gamma,P1,P0);
    %rhs2 = jumpMu/mu0 * (0.5*M_P1_P0 - K') * HJn_coeffs;
    rhs2 = 1/mu0*(K' - 0.5 * M_P1_P0) * alpha_coeffs;
    rhs = [rhs1;rhs2;0];

    sol= blockmat\rhs;

    psi = sol(1:NP0);
    g = sol(NP0+1:NP0+NP1);
    
end