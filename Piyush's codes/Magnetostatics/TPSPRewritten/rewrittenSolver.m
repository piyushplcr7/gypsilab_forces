function [Tnu,Tdu] = rewrittenSolver(Gamma,omega_src,J,mu,mu0)
    % BEM spaces
    bndmesh = Gamma.msh;
    P1 = fem(bndmesh,'P1');
    P0 = fem(bndmesh,'P0');
    V = single_layer(Gamma,P0,P0);
    K = double_layer_laplace(Gamma,P0,P1);
    W = single_layer(Gamma,P1.nxgrad,P1.nxgrad);
    % Vector to enforce zero P1 mean constraint
    Vec = integral(Gamma,P1);

    % Computing the vector potential Hj at Gamma
    [X_bndmesh,~] = Gamma.qud;
    [Y_src,Wsrc] = omega_src.qud;
    NX = size(X_bndmesh,1);
    NY = size(Y_src,1);
    YY = repmat(Y_src,NX,1);
    XX = repelem(X_bndmesh,NY,1);
    JYY = J(YY(:,1),YY(:,2),YY(:,3));
    gradxG = 1/4/pi * (YY-XX)./vecnorm(XX-YY,2,2).^3;
    integrand = cross(gradxG,JYY,2);
    HJ1 = sum(Wsrc.*reshape(integrand(:,1),NY,NX),1)';
    HJ2 = sum(Wsrc.*reshape(integrand(:,2),NY,NX),1)';
    HJ3 = sum(Wsrc.*reshape(integrand(:,3),NY,NX),1)';
    normals = Gamma.qudNrm;
    HJdotnatGamma = dot([HJ1 HJ2 HJ3],normals,2);
    HJdotncoeffs = proj(HJdotnatGamma,Gamma,P0);

    jmu = mu0-mu;

    %% Solving linear system
    NP1 = P1.ndof;
    NP0 = P0.ndof;
    % LHS matrix
    blockmat = [(1+mu0/mu)*V -2*K zeros(NP0,1);
                2*K' (1+mu/mu0)*W Vec;
                zeros(1,NP0) Vec' 0];

    % RHS
    M01 = mass_matrix(Gamma,P0,P1);
    rhs1 = -jmu/mu * V *HJdotncoeffs;
    rhs2 = jmu/2/mu0 * M01'*HJdotncoeffs -jmu/mu0 * K' * HJdotncoeffs;
    rhs = [rhs1;rhs2;0];

    sol= blockmat\rhs;

    Tnu = sol(1:NP0);
    Tdu = sol(NP0+1:NP0+NP1);
end