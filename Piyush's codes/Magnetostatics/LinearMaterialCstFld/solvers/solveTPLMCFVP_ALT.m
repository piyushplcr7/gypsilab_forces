% Solver linear material constant field vector potential
% Solve Transmission Problem. Returns the exterior traces 

function [Psi_i,g_i,Psi_e] = solveTPLMCFVP_ALT(bndmesh_i,bndmesh_e,mu,mu_0,g_e)
    %% BEM Spaces
    Gamma_i = dom(bndmesh_i,3);
    % [X_i,~] = Gamma_i.qud;

    Gamma_e = dom(bndmesh_e,3);
    % [X_e,~] = Gamma_e.qud;

    % BEM spaces
    % curl conforming -> Dirichlet trace
    % NED_i = fem(bndmesh_i,'NED'); 
    P1_i = fem(bndmesh_i,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0_i = nxgrad(P1_i); 
    % Kernel of the surface curl operator
    % Div conforming space 
    DIV_i = fem(bndmesh_i,'RWG');

    NED_e = fem(bndmesh_e,'NED');
    P1_e = fem(bndmesh_e,'P1');
    DIV0_e = nxgrad(P1_e);
    DIV_e = fem(bndmesh_e,'RWG');

    %% Linear System

    blockopr = LMCFVP_systemMat(bndmesh_i,bndmesh_e,mu,mu_0);

    %% Computing the RHS

    Cei = double_layer_magnetostatics_cross(Gamma_i,Gamma_e,DIV0_i,DIV_e);
    Nei = -single_layer_cross(Gamma_i,Gamma_e,DIV_i.div,DIV_e.div);
    Cee = double_layer_magnetostatics(Gamma_e,DIV0_e,DIV_e);

    % rhs1
    rhs1 = Cei * g_e;

    % rhs2
    rhs2 = -Nei * g_e;

    % rhs3
    rhs3 = Cee * g_e + 0.5 * mass_matrix(Gamma_e,DIV0_e,NED_e) * g_e;

    rhs = [rhs1; rhs2; rhs3 ;0; zeros(P1_i.ndof,1); 0; 0];
    % rhs = [rhs1; rhs2; rhs3];
    sol = blockopr\rhs;

    % tol = 1e-12;
    % maxit = 200;

    % sol = gmres(blockopr,rhs,[],tol,maxit);

    Psi_i = sol(1:P1_i.ndof);
    g_i = sol(P1_i.ndof+1:P1_i.ndof+DIV_i.ndof);
    Psi_e = sol(P1_i.ndof+DIV_i.ndof+1: P1_i.ndof+DIV_i.ndof + P1_e.ndof);
end