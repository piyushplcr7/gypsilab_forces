function val = dsl3_ALT(bndmesh_i,bndmesh_e,g_e,g_i,Vel,DVel)
    P1_e = fem(bndmesh_e,'P1');
    P1_i = fem(bndmesh_i,'P1');

    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    nxgrad_g_e_vals = reconstruct(g_e,Gamma_e,nxgrad(P1_e));
    nxgrad_g_i_vals = reconstruct(g_i,Gamma_i,nxgrad(P1_i));

    [X,WX] = Gamma_i.qud;
    [Y,WY] = Gamma_e.qud;

    NX = size(X,1); NY = size(Y,1);

    % XX -> repelem (NY) ; YY -> repmat (NX)
    XX = repelem(X,NY,1); YY = repmat(Y,NX,1);
    W = repelem(WX,NY,1).*repmat(WY,NX,1);

    dsl3kernel = 1/4/pi ./vecnorm(XX-YY,2,2);

    nxgrad_g_e_vals_YY = repmat(nxgrad_g_e_vals,NX,1);

    DVel1 = DVel{1}(X);
    DVel2 = DVel{2}(X);
    DVel3 = DVel{3}(X);

    DVel_nxgrad_g_i = [dot(DVel1,nxgrad_g_i_vals,2) dot(DVel2,nxgrad_g_i_vals,2) dot(DVel3,nxgrad_g_i_vals,2)];
    DVel_nxgrad_g_i_XX = repelem(DVel_nxgrad_g_i,NY,1);

    val = -sum(W.* dsl3kernel .* dot(nxgrad_g_e_vals_YY,DVel_nxgrad_g_i_XX,2),1);

end