% Magnetostatics BEM implementation test

addpath(genpath("../../"));
clear; clc;

ivals = 6:11;
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
    L = 5*[1 1 1];
    T = [0 0 0];
    
    % Solution domain
    bndmesh = bndmeshCubeTranslated(N,L,T);
    % Spherical mesh
    %bndmesh = mshSphere(N,1);
    hvals(i) = sqrt(mean(bndmesh.ndv,1));

    Gamma = dom(bndmesh,7);

    %% BEM Galerkin matrices
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

    %% Generating synthetic traces 
    % Evaluation points on Gamma
    [X,Wt,elt2qud] = Gamma.qud;

    % Taking traces of the fields
    normals = Gamma.qudNrm;
    Nnrm = size(normals,1);
    A = zeros(Nnrm,3);
    A(:,1) = 1+A(:,1);
    TDA = A - dot(A,normals,2).*normals;

    % Projecting the Dirichlet trace to NED space
    TDA_NED_coeffs = proj(TDA,Gamma,NED);

    %% Solving for the Neumann Trace
    % Creating the rhs vector
    % RHS
    rhs = (0.5*Mmat+Cmat)*TDA_NED_coeffs;
    rhs = [rhs;0];
    norm(rhs)
    sol = Amod\rhs;
    TNA_sol_DIV0_coeffs = sol(1:end-1);
    
    %% Computing the errors
    
    % Need to project computed Neumann trace to DIV space
    TNA_sol = reconstruct(TNA_sol_DIV0_coeffs,Gamma,DIV0);

    TNA_sol_DIV_coeffs = proj(TNA_sol,Gamma,DIV);

    MDD = mass_matrix(Gamma,DIV,DIV);
    err_DIV_coeffs = TNA_sol_DIV_coeffs;
    L2err = err_DIV_coeffs'*MDD*err_DIV_coeffs
    Hdiverr = err_DIV_coeffs'*single_layer(Gamma,DIV,DIV)*err_DIV_coeffs
    Hdiv0err = TNA_sol_DIV0_coeffs'*Amat*TNA_sol_DIV0_coeffs
    L2errs(i) = L2err;
    Hdiverrs(i) = Hdiv0err;

    if ~true
        quiver3wrapper(X,TNA,'blue');
        hold on;
        quiver3wrapper(X,TNA_sol,'red');
    end
end
loglog(hvals,Hdiverrs,'-s');
hold on;
loglog(hvals,L2errs,'-*')
%save('cleaned_constant_td.mat','L2errs','Hdiverrs','hvals');