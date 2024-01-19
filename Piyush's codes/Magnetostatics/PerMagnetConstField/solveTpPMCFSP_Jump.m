function [psi_i,g_i,psi_e] = solveTpPMCFSP_Jump(bndmesh_i,bndmesh_e,jump_coeffs)
    % BEM Spaces
    P1_i = fem(bndmesh_i,'P1');
    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');
    order = 7;
    Gamma_i = dom(bndmesh_i,order);
    Gamma_e = dom(bndmesh_e,order);

    %% Operator matrices

    Vii = single_layer(Gamma_i,P0_i,P0_i);
    Vee = single_layer(Gamma_e,P0_e,P0_e);

    Kii = double_layer_laplace(Gamma_i,P0_i,P1_i);
    Wii = single_layer(Gamma_i,P1_i.nxgrad,P1_i.nxgrad);
    
    % Cross matrices
    Vei = single_layer_cross(Gamma_i,Gamma_e,P0_i,P0_e);
    Kie = double_layer_laplace_cross(Gamma_e,Gamma_i,P0_e,P1_i);

    %% Linear system

    blockopr = [2*Vii, 2*Kii, Vei;
                -2*Kii', 2*Wii, -Kie';
                Vei', Kie, Vee];

    rhs = [-Vii * jump_coeffs;
            -0.5 * mass_matrix(Gamma_i,P1_i,P0_i) * jump_coeffs + Kii' * jump_coeffs;
            zeros(P0_e.ndof,1)];

    sol = blockopr\rhs;

    psi_i = sol(1:P0_i.ndof);
    g_i = sol(P0_i.ndof+1:P0_i.ndof+P1_i.ndof);
    psi_e = sol(P0_i.ndof+P1_i.ndof+1:P0_i.ndof+P1_i.ndof+P0_e.ndof);
end