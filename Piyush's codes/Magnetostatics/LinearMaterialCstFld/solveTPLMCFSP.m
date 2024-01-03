function [psi_i,g_i,psi_e,psi_I_recon,g_i_recon] = solveTPLMCFSP(bndmesh_i,bndmesh_e,mu,mu0,H0)
    % BEM Spaces
    P1_i = fem(bndmesh_i,'P1');
%     P1_e = fem(bndmesh_e,'P1');

    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');

    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    normals_i = Gamma_i.qudNrm;

    jumpMu = mu0-mu;

    % Computing H0 . n at Gamma_i
    H0dotn = normals_i * H0';
    H0dotncoeffs = proj(H0dotn,Gamma_i,P0_i);

    %% Operator matrices

    Vii = single_layer(Gamma_i,P0_i,P0_i);
    Vee = single_layer(Gamma_e,P0_e,P0_e);

    Kii = double_layer_laplace(Gamma_i,P0_i,P1_i);
    Wii = single_layer(Gamma_i,P1_i.nxgrad,P1_i.nxgrad);

    Vec = integral(Gamma_i,P1_i);
    
    % Cross matrices
    Vei = single_layer_cross(Gamma_i,Gamma_e,P0_i,P0_e);
    Kie = double_layer_laplace_cross(Gamma_e,Gamma_i,P0_e,P1_i);


    %% Linear system

    blockopr = [(1+mu0/mu)*Vii, 2*Kii, Vei, zeros(P0_i.ndof,1);
                -2*Kii', (1+mu/mu0)*Wii, -Kie', Vec;
                Vei', Kie, Vee, zeros(P0_e.ndof,1);
                zeros(1,P0_i.ndof), Vec', zeros(1,P0_e.ndof), 0];

    rhs = [jumpMu/mu * Vii * H0dotncoeffs;
            jumpMu/2/mu0 * mass_matrix(Gamma_i,P1_i,P0_i) * H0dotncoeffs - jumpMu/mu0 * Kii' * H0dotncoeffs;
            zeros(P0_e.ndof,1);
            0];

    sol = blockopr\rhs;

    psi_i = sol(1:P0_i.ndof);
    g_i = sol(P0_i.ndof+1:P0_i.ndof+P1_i.ndof);
    psi_e = sol(P0_i.ndof+P1_i.ndof+1:P0_i.ndof+P1_i.ndof+P0_e.ndof);

    %% Finding the traces of the solution 

    g_I_star = g_i;
    psi_I_star = jumpMu/mu * H0dotncoeffs - mu0/mu * psi_i;

    M01 = mass_matrix(Gamma_i,P0_i,P1_i);
    M10 = M01';
    
    g_i_recon = M01\(Vii * psi_I_star + 0.5 * M01 * g_I_star - Kii * g_I_star);

    psi_I_recon = M10\(0.5 * M10 * psi_I_star + Kii' * psi_I_star + Wii * g_I_star);

    
end