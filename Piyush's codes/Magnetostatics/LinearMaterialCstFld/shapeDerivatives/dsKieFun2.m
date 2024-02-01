% Outer integral (test) on Gamma_e <-> X. 
% Inner integral (trial) on Gamma_i <-> Y
function val = dsKieFun2(bndmesh_i,bndmesh_e,g_i,psi_e,Vel,DVel)
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

    DVel1iYY = DVel{1}(YY);
    DVel2iYY = DVel{2}(YY);
    DVel3iYY = DVel{3}(YY);
    divVeliYY = DVel1iYY(:,1) + DVel2iYY(:,2) + DVel3iYY(:,3);
    
    gradyG = 1/4/pi * (XX-YY)./vecnorm(XX-YY,2,2).^3;
    DVelgradyG = [dot(DVel1iYY,gradyG,2) dot(DVel2iYY,gradyG,2) dot(DVel3iYY,gradyG,2)];
    dskie2kernel = divVeliYY.*dot(gradyG,normalsi_YY,2) - dot(normalsi_YY,DVelgradyG,2);
    dsKie2 = sum(W.*dskie2kernel.*g_i_vals_YY.*psi_e_vals_XX,1);

    val = dsKie2;

end