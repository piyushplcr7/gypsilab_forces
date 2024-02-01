

function val = dsl1_ALT(bndmesh_i,bndmesh_e,g_e,psi_i,Vel,DVel)
    P1_e = fem(bndmesh_e,'P1');
    P0_i = fem(bndmesh_i,'P0');

    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    g_e_vals = reconstruct(g_e,Gamma_e,P1_e);
    psi_i_vals = reconstruct(psi_i,Gamma_i,P0_i);

    [X,WX] = Gamma_i.qud;
    [Y,WY] = Gamma_e.qud;

    NX = size(X,1); NY = size(Y,1);

    % XX -> repelem (NY) ; YY -> repmat (NX)
    XX = repelem(X,NY,1); YY = repmat(Y,NX,1);
    W = repelem(WX,NY,1).*repmat(WY,NX,1);

    normals_e = Gamma_e.qudNrm;

    normals_e_YY = repmat(normals_e,NX,1);

    dsl1kernel = 1/4/pi * (...
    3 * dot(XX-YY,Vel(XX),2).*dot(XX-YY,normals_e_YY,2)./vecnorm(XX-YY,2,2).^5 ...
    - dot(Vel(XX),normals_e_YY,2)./vecnorm(XX-YY,2,2).^3 );

    g_e_vals_YY = repmat(g_e_vals,NX,1);
    psi_i_vals_XX = repelem(psi_i_vals,NY,1);

    val = -sum(W.* dsl1kernel .* g_e_vals_YY.*psi_i_vals_XX,1);

end