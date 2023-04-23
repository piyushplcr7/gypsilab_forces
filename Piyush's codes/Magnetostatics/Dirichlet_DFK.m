function [TNA_sol_DIV0_coeffs] = Dirichlet_DFK(TDA,DIV0,order)
    bndmesh = DIV0.msh;
    Gamma = dom(bndmesh,order);

    % BEM spaces
    NED = fem(bndmesh,'NED'); % curl conforming -> Dirichlet trace
    P1 = fem(bndmesh,'P1');
    DIV0 = nxgrad(P1); % div conforming, div0 -> Neumann trace
    DIV = fem(bndmesh,'RWG');

    % For operator A
    Amat = single_layer(Gamma,DIV0,DIV0);
    % For operator C
    Cmat = double_layer_magnetostatics(Gamma,DIV0,DIV);
    Mmat = mass_matrix(Gamma,DIV0,NED);
    
    % Regularize by removing constants
    % Vector to enforce zero mean for functions
    B = integral(Gamma,P1);
    Amod = [Amat B; B' 0];

    % Projecting the Dirichlet trace to NED space
    TDA_NED_coeffs = proj(TDA,Gamma,NED);
    
    rhs = (0.5*Mmat+Cmat)*TDA_NED_coeffs;
    rhs = [rhs;0];
    sol = Amod\rhs;
    TNA_sol_DIV0_coeffs = sol(1:end-1);
end