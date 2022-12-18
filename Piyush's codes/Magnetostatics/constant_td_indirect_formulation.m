% Magnetostatics BEM implementation test

addpath(genpath("../../"));
clear; clc;

ivals = 6:6;
%ivals = 100:500:8600
Nivals = size(ivals,2);

hvals = ivals*0;

% Projection errors
TD_proj_errs = zeros(Nivals,1);
TN_proj_errs = TD_proj_errs;

% Solution errors
L2errs = TD_proj_errs;
Hdiverrs = L2errs*0;

for i = 1:Nivals
    N = 2^ivals(i)
    %N = ivals(i)
    % Cube size and position
    L = [1 1 1];
    
    volmesh = mshCube(N,L);
    % Solution domain
    bndmesh = volmesh.bnd;
    % Spherical mesh
    %bndmesh = mshSphere(N,1);
    hvals(i) = sqrt(mean(bndmesh.ndv,1));

    Gamma = dom(bndmesh,3);

    %% BEM Galerkin matrices
    % BEM spaces
    NED = fem(bndmesh,'NED'); % curl conforming -> Dirichlet trace
    P1 = fem(bndmesh,'P1');
    DIV0 = nxgrad(P1); % div conforming, div0 -> Neumann trace
    DIV = fem(bndmesh,'RWG');
    

    W = hypersingular_laplace(Gamma,P1,P1);
    Cmat = double_layer(Gamma,DIV0,NED);
    Mmat = mass_matrix(Gamma,DIV0,NED);
    
    % Regularize by removing constants
    % Vector to enforce zero mean for functions
    B = integral(Gamma,P1);
    Amod = [W B; B' 0];

    %% Generating synthetic traces 
    % Evaluation points on Gamma
    [X,Wt,elt2qud] = Gamma.qud;

    % Taking traces of the fields
    normals = Gamma.qudNrm;
    A = ones(size(normals));
    TDA = A - dot(A,normals,2).*normals;

    % Projecting the Dirichlet trace to NED space
    TDA_NED_coeffs = proj(TDA,Gamma,NED);

    %% Solving for the Ansatz density
    % Creating the rhs vector
    % RHS
    rhs = Mmat*TDA_NED_coeffs;
    rhs = [rhs;0];
    sol = Amod\rhs;
    sol_DIV0_coeffs = sol(1:end-1);
    
    %% Evaluating the ansatz at interior points to check solution
    Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); 
    Xi = volmesh.ctr;
    eval_A = integral(Xi,Gamma,Gxy,DIV0);

   
end
%save('cleaned_constant_td.mat','L2errs','Hdiverrs','hvals');