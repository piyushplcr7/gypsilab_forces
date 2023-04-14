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

    % BEM spaces
    NED = fem(bndmesh,'NED'); % curl conforming -> Dirichlet trace
    P1 = fem(bndmesh,'P1');
    DIV0 = nxgrad(P1); % div conforming, div0 -> Neumann trace
    DIV = fem(bndmesh,'RWG');

    Gamma = dom(bndmesh,3);

    %% Generating synthetic traces from a source
    %[TDA,TNA] = getTraces1(Gamma);
    [TDA,TNA] = getTraces2(Gamma);
    %[TDA,TNA,mesh_src] = getTracesTorusSource(Gamma);
    TNA_DIV_coeffs = proj(TNA,Gamma,DIV);

    %% Solving for the Neumann Trace
    TNA_sol_DIV0_coeffs = Dirichlet_DFK(TDA,DIV0);
    
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
    [X,~] = Gamma.qud;
    if true
        quiver3wrapper(X,TNA,'blue');
        hold on;
        quiver3wrapper(X,TNA_sol,'red');
    end
end
figure;
loglog(hvals,Hdiverrs,'-s');
hold on;
loglog(hvals,L2errs,'-*')

save('cleaned.mat','L2errs','Hdiverrs','hvals');