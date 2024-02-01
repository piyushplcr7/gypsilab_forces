

function val = dsl2_ALT(bndmesh_i,bndmesh_e,g_e,g_i,Vel,DVel)
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

    dsl2kernel = 1/4/pi * dot(YY-XX,Vel(XX),2)./vecnorm(XX-YY,2,2).^3;

    nxgrad_g_e_vals_YY = repmat(nxgrad_g_e_vals,NX,1);
    nxgrad_g_i_vals_XX = repelem(nxgrad_g_i_vals,NY,1);

    val = -sum(W.* dsl2kernel .* dot(nxgrad_g_e_vals_YY,nxgrad_g_i_vals_XX,2),1);

end