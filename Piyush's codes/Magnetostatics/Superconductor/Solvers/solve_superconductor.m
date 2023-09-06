% Solver for superconducting domain. Zero Dirichlet boundary conditions 
% Imposed

% Inputs
% bndmesh : Boundary mesh for the superconducting domain
% J : Function object that takes in vectorized inputs (x,y,z) and returns
%     the source current in NX3 format
% mesh_src : Mesh for the current source

% Output
% TnA0 : Neumann trace in div0 space
% TnA : Neumann trace in RWG space

function [TnA0,TnA] = solve_superconductor(bndmesh, J, mesh_src)
    % BEM Spaces
    % curl conforming -> Dirichlet trace
    NED = fem(bndmesh,'NED'); 
    P1 = fem(bndmesh,'P1');
    % Div conforming space 
    DIV = fem(bndmesh,'RWG');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0 = nxgrad(P1); 

    % Dom objects for integration
    Gamma = dom(bndmesh,3);
    omega_src = dom(mesh_src,3);

    % Galerkin matrix for A
    Amat = single_layer(Gamma,DIV0,DIV0);

    %% Computing the RHS

    % Evaluation points for the Dirichlet Trace
    [X,~,~] = Gamma.qud;

    % Computing the Newton Potential on the points X
    Asrc = compute_vecpot(J,omega_src,X);
    % Taking traces of the fields at the points X
    normals = Gamma.qudNrm;
    TDAsrc = Asrc - dot(Asrc,normals,2).*normals;

    % Projecting the Dirichlet trace to NED space
    TDAsrc_NED_coeffs = proj(TDAsrc,Gamma,NED);

    % RHS
    rhs = mass_matrix(Gamma,DIV0,NED) * TDAsrc_NED_coeffs;

    %% Blockmat for LHS
    % Vector to enforce zero mean for P1 functions
    vec = integral(Gamma,P1);
    blockmat = [Amat vec;
                vec' 0];
    rhsmodif = [rhs;0];

    %% Solving the linear system
    solextended = blockmat\rhsmodif;
    TnA0 = solextended(1:end-1);

    TnA_sol = reconstruct(TnA0,Gamma,DIV0);
    TnA = proj(TnA_sol,Gamma,DIV);

end

