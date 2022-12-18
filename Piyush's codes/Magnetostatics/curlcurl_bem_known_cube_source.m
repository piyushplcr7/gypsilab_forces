% Magnetostatics BEM implementation test

addpath(genpath("../../"));

clear; clc;

ivals = 6:6;
Nivals = size(ivals,2);
TD_proj_errs = zeros(Nivals,1);
TN_proj_errs = TD_proj_errs;
L2errs = TD_proj_errs;
Hdiverrs = L2errs*0;

plot = ~true;

for i = 1:Nivals
    N = 2^ivals(i)
    
    % Cube size and position
    L = [1 1 1];
    T = [5 0 0];
    
    % Solution domain
    bndmesh = bndmeshCubeTranslated(N,L,T);

    Gamma = dom(bndmesh,3);

    %% BEM Galerkin matrices

    % BEM spaces
    NED = fem(bndmesh,'NED'); % curl conforming -> Dirichlet trace
    P1 = fem(bndmesh,'P1');
    DIV0 = nxgrad(P1); % div conforming, div0 -> Neumann trace
    DIV = fem(bndmesh,'RWG');

    Amat = single_layer(Gamma,DIV0,DIV0);
    Cmat = double_layer(Gamma,DIV0,NED);
    Mmat = mass_matrix(Gamma,DIV0,NED);

    % The matrix A on LHS is singular, regularize it
    alpha = 10;
    Amod = Amat + alpha * hypersingular_regularizer(bndmesh,3);

    %% Generating synthetic traces from a source
    N_src = N;
    mesh_src = mshCube(N_src,[1 1 1]);
    mesh_src = mesh_src.translate([0.5 0.5 0.5]);
    omega_src = dom(mesh_src,3);

    % Vector potential
    Aknown = @(x,y,z) [cos(pi*x).*sin(pi*y).*sin(pi*z),...
             -0.5*sin(pi*x).*cos(pi*y).*sin(pi*z),...
             -0.5*sin(pi*x).*sin(pi*y).*cos(pi*z)];

    % Source 
    J = @(x,y,z) -3*pi^2*Aknown(x,y,z);

    % Evaluation points on Gamma
    [X,W,elt2qud] = Gamma.qud;

    % Computing the fields on the points X
    A = compute_vecpot(J,omega_src,X);
    curlA = compute_vecpot_curl(J,omega_src,X);

    % Taking traces of the fields
    normals = Gamma.qudNrm;
    TDA = A - dot(A,normals,2).*normals;
    TNA = cross(normals,curlA,2);

    % Projecting the Dirichlet trace to NED space
    TDA_NED_coeffs = proj(TDA,Gamma,NED);
    TDA_NED_recon = reconstruct(TDA_NED_coeffs,Gamma,NED);
    L2err_TDA_proj = sum(W.*(vecnorm(TDA-TDA_NED_recon,2,2).^2),1)

    TD_proj_errs(i) = L2err_TDA_proj;

    % Projecting the Neumann trace to DIV space
    TNA_DIV_coeffs = proj(TNA,Gamma,DIV);
    TNA_DIV_recon = reconstruct(TNA_DIV_coeffs,Gamma,DIV);
    L2err_TNA_proj = sum(W.*(vecnorm(TNA-TNA_DIV_recon,2,2).^2),1)

    TN_proj_errs(i) = L2err_TNA_proj;

    %% Solving for the Neumann Trace
    TNA_sol_DIV0_coeffs = Amod\((0.5*Mmat+Cmat)*TDA_NED_coeffs);
    
    %% Computing the errors
    
    % Need to project computed Neumann trace to DIV space
    TNA_sol = reconstruct(TNA_sol_DIV0_coeffs,Gamma,DIV0);

    TNA_sol_DIV_coeffs = proj(TNA_sol,Gamma,DIV);

    MDD = mass_matrix(Gamma,DIV,DIV);
    err_DIV_coeffs = TNA_DIV_coeffs-TNA_sol_DIV_coeffs;
    L2err = err_DIV_coeffs'*MDD*err_DIV_coeffs
    Hdiverr = err_DIV_coeffs'*single_layer(Gamma,DIV,DIV)*err_DIV_coeffs

    L2errs(i) = L2err;
    Hdiverrs(i) = Hdiverr;

    if plot
        quiver3wrapper(X,TNA,'blue');
        hold on;
        quiver3wrapper(X,TNA_sol,'red');
    end
end



