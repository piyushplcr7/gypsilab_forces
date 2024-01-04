function sd = sdBEMTpDielSC_ConstVel(bndmesh,epsilon,epsilon0,Tnu,Tdu,rho,omega_src,Vel,DVel)
    P0 = fem(bndmesh,'P0');
    P1 = fem(bndmesh,'P1');

    Gamma = dom(bndmesh,3);
    normals = Gamma.qudNrm;
    normals_src = omega_src.qudNrm;

    psi_vals = reconstruct(Tnu,Gamma,P0);
    g_vals = reconstruct(Tdu,Gamma,P1);

    [Y,WY] = omega_src.qud;
    [X,WX] = Gamma.qud;
    NY = size(Y,1);
    NX = size(X,1);

    YY = repmat(Y,NX,1);
    XX = repelem(X,NY,1);

    W = repmat(WY,NX,1).*repelem(WX,NY,1);

    gradxGdotvelx = 1/4/pi * dot(YY-XX,Vel(XX),2)./vecnorm(XX-YY,2,2).^3;
    rhoYY = rho(YY);
    psiXX = repelem(psi_vals,NY,1);
    gXX = repelem(g_vals,NY,1);

    t1 = - sum(W.* gradxGdotvelx.*rhoYY.*psiXX ,1);

    normals_src_YY = repmat(normals_src,NX,1);

    t2kernel = 3/4/pi * dot(XX-YY,normals_src_YY,2).*dot(XX-YY,Vel(XX),2)./vecnorm(XX-YY,2,2).^5 - 1/4/pi * dot(Vel(XX),normals_src_YY,2)./vecnorm(XX-YY,2,2).^3;
    t2 = sum(W.*t2kernel.*rhoYY.*gXX,1);

    sd = t1 + t2;


end