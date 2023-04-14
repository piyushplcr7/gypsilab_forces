% Magnetostatics BEM implementation test

addpath(genpath("../../"));
clear; clc; close all;

ivals = 6:6;
Nivals = size(ivals,2);

hvals = ivals*0;

% Projection errors
TD_proj_errs = zeros(Nivals,1);
TN_proj_errs = TD_proj_errs;

% Solution errors
L2errs = TD_proj_errs;
Hdiverrs = L2errs*0;

%plot = true;

for i = 1:Nivals
    N = 2^ivals(i)
    
    % Cube size and position
    L = 2*[1 1 1];
    T = [5 5 3];
    
    % Solution domain
    %bndmesh = bndmeshCubeTranslated(N,L,T);
    % Spherical mesh
    bndmesh = mshSphere(N,1);
    bndmesh = bndmesh.translate(T);

    hvals(i) = sqrt(mean(bndmesh.ndv,1));

    Gamma = dom(bndmesh,3);

    %% BEM Galerkin matrices

    % BEM spaces
    NED = fem(bndmesh,'NED'); % curl conforming -> Dirichlet trace
    P1 = fem(bndmesh,'P1');
    DIV0 = nxgrad(P1); % div conforming, div0 -> Neumann trace
    DIV = fem(bndmesh,'RWG');

    % For operator N
    Nmat = -single_layer(Gamma,DIV.div,DIV.div);
    % For operator B we can use C'
    Cmat0 = double_layer_magnetostatics(Gamma,DIV0,DIV);
    Cmat = double_layer_magnetostatics(Gamma,DIV,DIV);
    Mmat0 = mass_matrix(Gamma,DIV0,NED);
    Mmat = mass_matrix(Gamma,DIV,NED);
    
    % Regularize by removing constants
    % Vector to enforce zero mean for functions
    %B = integral(Gamma,P1);
    %Amod = [Amat B; B' 0];
    Amod = Nmat;

    %% Generating synthetic traces from a source
    N_src = N;
    R0 = 2;
    r0 = .5;
    [J,mesh_src] = get_torus_source(N_src,R0,r0);
    omega_src = dom(mesh_src,3);

    % Evaluation points on Gamma
    [X,Wt,elt2qud] = Gamma.qud;

    % Computing the fields on the points X
    A = compute_vecpot(J,omega_src,X);
    curlA = compute_vecpot_curl(J,omega_src,X);

    % Taking traces of the fields
    normals = Gamma.qudNrm;
    TDA = A - dot(A,normals,2).*normals;
    TNA = cross(curlA,normals,2);

    % Projecting the Dirichlet trace to NED space
    TDA_NED_coeffs = proj(TDA,Gamma,NED);

    % Projecting the Neumann trace to DIV space
    TNA_DIV_coeffs = proj(TNA,Gamma,DIV);
    TNA_DIV0_coeffs = proj(TNA,Gamma,DIV0);

    %% Solving for the Neumann Trace
    % Creating the rhs vector
    % RHS
    rhs = (-0.5*Mmat'+Cmat')*TNA_DIV_coeffs;
    %rhs = [rhs;0];
    sol = Amod\rhs;
    %TNA_sol_DIV0_coeffs = sol(1:end-1);
    TDA_sol_DIV_coeffs = sol;
    
    %% Computing the errors
    
    % Need to project computed Neumann trace to DIV space
    TDA_sol = reconstruct(TDA_sol_DIV_coeffs,Gamma,DIV);

    if true
        quiver3wrapper(X,TDA,'blue');
        hold on;
        quiver3wrapper(X,TDA_sol,'red');
    end
end
