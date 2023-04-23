% Solving for unknown Neumann trace given a Dirichlet trace using the 
% Second kind BIE

% TDA is passed as a vector of values of the neumann trace at qud pts of 
% Gamma

% Neumann solution sought in RWG space to make a square matrix

function [TNA] = Dirichlet_DSK(TDA,RWG,order)
    bndmesh = RWG.msh;
    Gamma = dom(bndmesh,order);

    % Curl conforming space
    NED = fem(bndmesh,'NED');

    %% LHS
    % Galerkin matrices

    % For operator N
    Nmat = -single_layer(Gamma,RWG.div,RWG.div);
  
    % For operator B, we can simply use C'
    % Operator C
    Cmat = double_layer_magnetostatics(Gamma,RWG,RWG);
    M = mass_matrix(Gamma,NED,RWG);
    TDA_coeffs = proj(TDA,Gamma,NED);

    lhsmat = (0.5*M - Cmat');
    rhs = -Nmat * TDA_coeffs;

    TNA = lhsmat\rhs;


end