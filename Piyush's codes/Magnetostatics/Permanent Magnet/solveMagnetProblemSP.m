% Solving transmission problem for permanent magnet using scalar potential

function [psi,g] = solveMagnetProblemSP(mesh,Gamma,mu,mu0,J,omega_src,M)
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

    % Computing the interacting bnd-vol integrals
    KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
    KDL = @(x,y,z)1/(4*pi)*z./vecnorm(z,2,2).^3;
    Nelt_bnd = bndmesh.nelt;
    Nelt_vol = mesh.nelt;
    RWGVol = fem(mesh,'RWG');
    omega = dom(mesh,4);
    [Xvol,~] = omega.qud;
    Mvol = M(Xvol);
    Mvol_coeffs = proj(Mvol,mesh,RWGVol);
    [II,JJ] = meshgrid(1:Nelt_bnd,1:Nelt_vol);
    intmat1 = panel_assembly_cross_modif(bndmesh,mesh,KV,RWGVol.div,P0,II(:),JJ(:));
    intmat2 = panel_assembly_cross_modif(bndmesh,mesh,KDL,RWGVol.div,P1,II(:),JJ(:));

    % RHS
    rhs1 = 1/mu * V * alpha_coeffs - intmat1 * Mvol_coeffs;
    M_P1_P0 = mass_matrix(Gamma,P1,P0);
    %rhs2 = jumpMu/mu0 * (0.5*M_P1_P0 - K') * HJn_coeffs;
    rhs2 = 1/mu0*(K' - 0.5 * M_P1_P0) * alpha_coeffs - mu/mu0 * intmat2 * Mvol_coeffs;
    rhs = [rhs1;rhs2;0];

    sol= blockmat\rhs;

    psi = sol(1:NP0);
    g = sol(NP0+1:NP0+NP1);
    
end