% Outer integral (test) on Gamma_e <-> X. 
% Inner integral (trial) on Gamma_i <-> Y
function val = dsKieFun1(bndmesh_i,bndmesh_e,g_i,psi_e,Vel,DVel)
    P1_i = fem(bndmesh_i,'P1');
    P0_e = fem(bndmesh_e,'P0');
    
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    g_i_vals = reconstruct(g_i,Gamma_i,P1_i);
    psi_e_vals = reconstruct(psi_e,Gamma_e,P0_e);

    
    [Y,WY] = Gamma_i.qud;
    [X,WX] = Gamma_e.qud;

    NX = size(X,1); NY = size(Y,1);
    % XX -> repelem (NY) ; YY -> repmat (NX)
    XX = repelem(X,NY,1); YY = repmat(Y,NX,1);
    W = repelem(WX,NY,1).*repmat(WY,NX,1);
    
    normals_i = Gamma_i.qudNrm;
    normalsi_YY = repmat(normals_i,NX,1);

    g_i_vals_YY = repmat(g_i_vals,NX,1);
    psi_e_vals_XX = repelem(psi_e_vals,NY,1);

    dskie1kernel = 1/4/pi * (3 * dot(XX-YY,normalsi_YY,2).*dot(XX-YY,Vel(YY),2)./vecnorm(XX-YY,2,2).^5 - dot(Vel(YY),normalsi_YY,2)./vecnorm(XX-YY,2,2).^3);
    dsKie1 = sum(W.*dskie1kernel.*g_i_vals_YY.*psi_e_vals_XX,1);


    val = dsKie1;

end