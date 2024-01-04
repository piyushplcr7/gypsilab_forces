function [Tnu,Tdu] = solveTpDielSourceCharge(bndmesh,epsilon,epsilon0,rho,omega_src)
    % BEM spaces
    P1 = fem(bndmesh,'P1');
    P0 = fem(bndmesh,'P0');

    Gamma = dom(bndmesh,3);

    V = single_layer(Gamma,P0,P0);
    K = double_layer_laplace(Gamma,P0,P1);
    W = single_layer(Gamma,P1.nxgrad,P1.nxgrad);
    % Vector to enforce zero P1 mean constraint
    Vec = integral(Gamma,P1);
    
    % Computing the potential due to source charge at points of bndmesh
    [X_bndmesh,~] = Gamma.qud;
    [Y_src,Wsrc] = omega_src.qud;
    normals_X_bndmesh = Gamma.qudNrm;
    NX = size(X_bndmesh,1);
    NY = size(Y_src,1);
    YY = repmat(Y_src,NX,1);
    XX = repelem(X_bndmesh,NY,1);
    G = 1/4/pi./vecnorm(XX-YY,2,2);
    G = reshape(G,NY,NX);
    gradxG = 1/4/pi * (YY-XX)./vecnorm(XX-YY,2,2).^3;
    normals_XX = repelem(normals_X_bndmesh,NY,1);
    gradxGdotnormalsXX = dot(gradxG,normals_XX,2);
    gradxGdotnormalsXX = reshape(gradxGdotnormalsXX,NY,NX);
    rho = rho(Y_src);
    
    TdNrho = sum(Wsrc.* rho .* G,1)';
    TnNrho = sum(Wsrc.* rho.* gradxGdotnormalsXX,1)';

    % Projecting these computed traces to appropriate spaces
    TdNrho_coeffs = proj(TdNrho,Gamma,P1);
    TnNrho_coeffs = proj(TnNrho,Gamma,P0);

    %% Solving linear system
    NP1 = P1.ndof;
    NP0 = P0.ndof;
    % LHS matrix
    blockmat = [(1+epsilon0/epsilon)*V -2*K zeros(NP0,1);
                2*K' (1+epsilon/epsilon0)*W Vec;
                zeros(1,NP0) Vec' 0];

    % RHS
    M01 = mass_matrix(Gamma,P0,P1);
    rhs1 = 1/epsilon0 * M01 * TdNrho_coeffs;
    rhs2 = 1/epsilon0 * M01' * TnNrho_coeffs;
    rhs = [rhs1;rhs2;0];

    sol= blockmat\rhs;

    Tnu = sol(1:NP0);
    Tdu = sol(NP0+1:NP0+NP1);


end