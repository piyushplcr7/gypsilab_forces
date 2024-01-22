function [psi_i,g_i,psi_e] = solveTPLMCFSP_ALT(bndmesh_i,bndmesh_e,mu,mu0,H0)
    % BEM Spaces
    P1_i = fem(bndmesh_i,'P1');
    P1_e = fem(bndmesh_e,'P1');

    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');

    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    normals_i = Gamma_i.qudNrm;

    jumpMu = mu0-mu;
    
    % Dirichlet trace at the outer boundary
    [X_e,~] = Gamma_e.qud;
    g_e_vals = X_e * H0';
    g_e_coeffs = proj(g_e_vals,Gamma_e,P1_e);

    % Computing H0 . n at Gamma_i
    H0dotn = normals_i * H0';
    H0dotncoeffs = proj(H0dotn,Gamma_i,P0_i);

    %% Operator matrices

    Vii = single_layer(Gamma_i,P0_i,P0_i);
    Vee = single_layer(Gamma_e,P0_e,P0_e);

    Kii = double_layer_laplace(Gamma_i,P0_i,P1_i);
    Kee = double_layer_laplace(Gamma_e,P0_e,P1_e);
    Wii = single_layer(Gamma_i,P1_i.nxgrad,P1_i.nxgrad);

    Vec = integral(Gamma_i,P1_i);
    
    % Cross matrices
    Vei = single_layer_cross(Gamma_i,Gamma_e,P0_i,P0_e);
    Kie = double_layer_laplace_cross(Gamma_e,Gamma_i,P0_e,P1_i);

    Kei = double_layer_laplace_cross(Gamma_i,Gamma_e,P0_i,P1_e);
    Wei = single_layer_cross(Gamma_i,Gamma_e,nxgrad(P1_i),nxgrad(P1_e));
    
    % Mass
    M01e = mass_matrix(Gamma_e,P0_e,P1_e);

    %% Linear system

    % blockopr = [(1+mu0/mu)*Vii, 2*Kii, Vei, zeros(P0_i.ndof,1);
    %             -2*Kii', (1+mu/mu0)*Wii, -Kie', Vec;
    %             Vei', Kie, Vee, zeros(P0_e.ndof,1);
    %             zeros(1,P0_i.ndof), Vec', zeros(1,P0_e.ndof), 0];

    blockopr = [(1+mu0/mu)*Vii, 2*Kii, Vei;
                -2*Kii', (1+mu/mu0)*Wii, -Kie';
                Vei', Kie, Vee];
    % cond(blockopr)

%     rhs = [Kei * g_e_coeffs;
%            Wei * g_e_coeffs;
%            Kee * g_e_coeffs + 0.5 * M01e * g_e_coeffs;
%             0];

    rhs = [Kei * g_e_coeffs;
           Wei * g_e_coeffs;
           Kee * g_e_coeffs + 0.5 * M01e * g_e_coeffs];

    sol = blockopr\rhs;

    psi_i = sol(1:P0_i.ndof);
    g_i = sol(P0_i.ndof+1:P0_i.ndof+P1_i.ndof);
    psi_e = sol(P0_i.ndof+P1_i.ndof+1:P0_i.ndof+P1_i.ndof+P0_e.ndof);

    %% Testing reconstructed traces of Interior Rep. Formula at interface
%     M01i = mass_matrix(Gamma_i,P0_i,P1_i);
%     tested_g_I = -mu0/mu * Vii * psi_i + 0.5 * M01i * g_i - Kii * g_i;
%     tested_psi_I = -mu0/2/mu * M01i' * psi_i - mu0/mu * Kii' * psi_i + Wii * g_i;
% 
%     % Exact traces for mu = mu0
%     expl_psi_i_vals = -normals_i * H0';
%     psi_i_exact = proj(expl_psi_i_vals,Gamma_i,P0_i);
%     [X_i,~] = Gamma_i.qud;
%     expl_g_i_vals = X_i * H0';
%     g_i_exact = proj(expl_g_i_vals,Gamma_i,P1_i);
% 
%     % Comparing
%     tested_g_i_exact = M01i * g_i_exact;
%     tested_psi_I_exact = -mu0/mu * M01i' * psi_i_exact;
% 
%     fprintf("Err_g, max = %f, norm = %f \n",max(abs(tested_g_i_exact-tested_g_I)),norm(tested_g_i_exact-tested_g_I)/size(tested_g_I,1));
%     fprintf("Err_psi, max = %f, norm = %f \n",max(abs(tested_psi_I_exact-tested_psi_I)),norm(tested_psi_I_exact-tested_psi_I)/size(tested_psi_I,1));

    % disp(sol(end));
end