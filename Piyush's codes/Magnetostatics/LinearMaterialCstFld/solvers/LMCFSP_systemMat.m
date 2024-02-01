function blockopr = LMCFSP_systemMat(bndmesh_i,bndmesh_e,mu,mu0)

    P1_i = fem(bndmesh_i,'P1');

    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');

    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);
    
    Vii = single_layer(Gamma_i,P0_i,P0_i);
    Vee = single_layer(Gamma_e,P0_e,P0_e);

    Kii = double_layer_laplace(Gamma_i,P0_i,P1_i);
    Wii = single_layer(Gamma_i,P1_i.nxgrad,P1_i.nxgrad);

    % Cross matrices
    Vei = single_layer_cross(Gamma_i,Gamma_e,P0_i,P0_e);
    Kie = double_layer_laplace_cross(Gamma_e,Gamma_i,P0_e,P1_i);
    
    blockopr = [(1+mu0/mu)*Vii, 2*Kii, Vei;
                -2*Kii', (1+mu/mu0)*Wii, -Kie';
                Vei', Kie, Vee];

end