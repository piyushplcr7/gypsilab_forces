% Solver linear material constant field vector potential
% Solve Transmission Problem. Returns the exterior traces 

function [Psi_i,g_i,Psi_e] = solveTPConstantFld(bndmesh_i,bndmesh_e,mu,mu_0,B_0)
    %% BEM Spaces
    Gamma_i = dom(bndmesh_i,3);
    [X_i,~] = Gamma_i.qud;

    Gamma_e = dom(bndmesh_e,3);
    [X_e,~] = Gamma_e.qud;

    % BEM spaces
    % curl conforming -> Dirichlet trace
    NED_i = fem(bndmesh_i,'NED'); 
    P1_i = fem(bndmesh_i,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0_i = nxgrad(P1_i); 
    % Kernel of the surface curl operator
    Ker_curl_i = grad(P1_i); 
    % Div conforming space 
    DIV_i = fem(bndmesh_i,'RWG');

    NED_e = fem(bndmesh_e,'NED');
    P1_e = fem(bndmesh_e,'P1');
    DIV0_e = nxgrad(P1_e);
    DIV_e = fem(bndmesh_e,'RWG');

    %% Galerkin matrices 
    % For operator A
    Aii = single_layer(Gamma_i,DIV0_i,DIV0_i);
    % For operator C
    Cii = double_layer_magnetostatics(Gamma_i,DIV0_i,DIV_i);
    % For operator B, we can simply use C'
    Bii = Cii';
    % For operator N
    Nii = -single_layer(Gamma_i,DIV_i.div,DIV_i.div);

    % vector to enforce zero mean for P1 functions
    veci = integral(Gamma_i,P1_i);

    ortho = single_layer(Gamma_i,NED_i,Ker_curl_i); % Uses my implementation.

    % Cross matrices
    Aei = single_layer_cross(Gamma_i,Gamma_e,DIV0_i,DIV0_e);
    Cie = double_layer_magnetostatics_cross(Gamma_e,Gamma_i,DIV0_e,DIV_i);

    Aee = single_layer(Gamma_e,DIV0_e,DIV0_e);
    vece = integral(Gamma_e,P1_e);

    %% Linear System

    %                 P1_i.ndof,          DIV_i.ndof,                  P1_e.ndof,                   1,                  P1_i.ndof,                   1,                   1; 
    blockopr = [(1+mu/mu_0)*Aii,               2*Cii,                        Aei,                veci, zeros(P1_i.ndof,P1_i.ndof),  zeros(P1_i.ndof,1),  zeros(P1_i.ndof,1); %P1_i.ndof
                         -2*Bii,    -(1+mu_0/mu)*Nii,                       Cie', zeros(DIV_i.ndof,1),                      ortho, zeros(DIV_i.ndof,1), zeros(DIV_i.ndof,1); % DIV_i.ndof
                           Aei',                 Cie,                        Aee,  zeros(P1_e.ndof,1), zeros(P1_e.ndof,P1_i.ndof),                vece,  zeros(P1_e.ndof,1); % P1_e.ndof
                          veci', zeros(1,DIV_i.ndof),         zeros(1,P1_e.ndof),                   0,         zeros(1,P1_i.ndof),                   0,                   0; % 1
     zeros(P1_i.ndof,P1_i.ndof),              ortho', zeros(P1_i.ndof,P1_e.ndof),  zeros(P1_i.ndof,1), zeros(P1_i.ndof,P1_i.ndof),  zeros(P1_i.ndof,1),                veci; % P1_i.ndof
             zeros(1,P1_i.ndof), zeros(1,DIV_i.ndof),                      vece',                   0,         zeros(1,P1_i.ndof),                   0,                   0; %1
             zeros(1,P1_i.ndof), zeros(1,DIV_i.ndof),         zeros(1,P1_e.ndof),                   0,                      veci',                   0,                   0; %1
             ];

    %% Computing the RHS

    % Computing B_0 X n at quadrature points on Gamma_i
    normals_i = Gamma_i.qudNrm;

    B_0xn = cross(repmat([B_0(1) B_0(2) B_0(3)],size(normals_i,1),1),normals_i,2);
    B_0xn_coeffs = proj(B_0xn,Gamma_i,DIV0_i);
    % Jump
    jumpMuInv = 1/mu_0 - 1/mu;

    % rhs1
    rhs1 = mu * jumpMuInv * Aii * B_0xn_coeffs;

    % rhs2
    massmat = mass_matrix(Gamma_i,NED_i,DIV0_i);
    rhs2 = -mu_0 * jumpMuInv * Cii' * B_0xn_coeffs + mu_0 * jumpMuInv/2 * massmat * B_0xn_coeffs;

    rhs = [rhs1; rhs2; zeros(P1_e.ndof,1) ;0; zeros(P1_i.ndof,1); 0; 0];
    sol = blockopr\rhs;

    Psi_i = sol(1:P1_i.ndof);
    g_i = sol(P1_i.ndof+1:P1_i.ndof+DIV_i.ndof);
    Psi_e = sol(P1_i.ndof+DIV_i.ndof+1: P1_i.ndof+DIV_i.ndof + P1_e.ndof);

    %% Verification of the solution via BIEs

    % (5.2.4)
%     err1 = -Aii * Psi_i - Aei * Psi_e - Cii * g_i + 0.5 * mass_matrix(Gamma_i,DIV0_i,NED_i) * g_i;
    err1 = (1+mu/mu_0) * Aii * Psi_i + Aei * Psi_e +2 *  Cii * g_i - mu * jumpMuInv * Aii * B_0xn_coeffs;


    % (5.2.5)
%     err2 = Bii * Psi_i + Cie' * Psi_e + Nii * g_i + 0.5 * mass_matrix(Gamma_i,NED_i,DIV0_i) * Psi_i; 
    err2 = -2*Bii * Psi_i -(1+mu_0/mu)* Nii * g_i - Cie' * Psi_e + mu_0 * jumpMuInv * Bii * B_0xn_coeffs - 0.5 * mu_0 * jumpMuInv * massmat * B_0xn_coeffs; 

    % (5.2.6)
    err3 = -Aei' * Psi_i - Aee * Psi_e - Cie * g_i;

    % (5.2.7) 
    %err4 = -C
end