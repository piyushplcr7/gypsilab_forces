% Outer integral (test) on Gamma_e <-> X. 
% Inner integral (trial) on Gamma_i <-> Y
function val = dsVieFun(bndmesh_i,bndmesh_e,psi_i,psi_e,Vel,DVel)
    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');
    
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    psi_i_vals = reconstruct(psi_i,Gamma_i,P0_i);
    psi_e_vals = reconstruct(psi_e,Gamma_e,P0_e);

    
    [Y,WY] = Gamma_i.qud;
    [X,WX] = Gamma_e.qud;

    NX = size(X,1); NY = size(Y,1);
    % XX -> repelem (NY) ; YY -> repmat (NX)
    XX = repelem(X,NY,1); YY = repmat(Y,NX,1);
    W = repelem(WX,NY,1).*repmat(WY,NX,1);

    psi_i_vals_YY = repmat(psi_i_vals,NX,1);
    psi_e_vals_XX = repelem(psi_e_vals,NY,1);


    dsviekernel = 1/4/pi * dot(XX-YY,Vel(YY),2)./vecnorm(XX-YY,2,2).^3;

    val = sum(W.*dsviekernel.*psi_i_vals_YY.*psi_e_vals_XX,1);

end