function sd = SdBemLMCFVP_ConstVel(bndmesh_i,bndmesh_e,Psi_i,g_i,Psi_e,Vel,DVel,mu0,mu,B0)
    jumpMuInv = 1/mu0-1/mu;
    % Integration domain
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    % BEM Spaces for Gamma_i
    P1_i = fem(bndmesh_i,'P1');
    DIV0_i = nxgrad(P1_i);
    RWG_i = fem(bndmesh_i,'RWG');

    % BEM Spaces for gamma_e
    P1_e = fem(bndmesh_e,'P1');
    DIV0_e = nxgrad(P1_e);
    RWG_e = fem(bndmesh_e,'RWG');

    Psi_i_vals = reconstruct(Psi_i,Gamma_i,RWG_i);
    Psi_e_vals = reconstruct(Psi_e,Gamma_e,RWG_e);
    nxgvals_i = reconstruct(g_i,Gamma_i,RWG_i);

    % Projecting B0xn to RWG_i space
    normals_i = Gamma_i.qudNrm;
    % B0xn at quadrature points
    B0xn = cross(repmat([B0(1) B0(2) B0(3)], size(normals_i,1),1),normals_i,2);
    B0xn_coeffs = proj(B0xn,Gamma_i,RWG_i);


    %% Non SS Computations 
    % Cross bilinear forms
    gradxGdotVelx = @(x,y) 1/4/pi ./vecnorm(x-y,2,2).^3 .* dot(y-x,Vel(x),2);

    % ei partial derivative computation
    Aei1mat = integral(Gamma_i,Gamma_e,RWG_i,gradxGdotVelx,RWG_e);
    % int_{Gamma_i} int_{\Gamma_e} gradxG(x,y).vel(x) psi_e(y).zeta_i(x) 
    Aei1 = Psi_i' * Aei1mat * Psi_e;

    % Cie computation
    % 1st term
    [X_e,WX_e] = Gamma_e.qud;
    [Y_i,WY_i] = Gamma_i.qud;
    NX = size(X_e,1);
    NY = size(Y_i,1);
    XX = repmat(X_e,NY,1); WWX = repmat(WX_e,NY,1); 
    YY = repelem(Y_i,NX,1); WWY = repelem(WY_i,NX,1);
    W = WWX .* WWY;

    nxgvals_iYY = repelem(nxgvals_i,NX,1);
    Psi_eXX = repmat(Psi_e_vals,NY,1);
    % Ciekernel1 = 3/4/pi * (XX-YY).*dot(XX-YY,Vel(XX)-Vel(YY),2)./vecnorm(XX-YY,2,2).^5 ...
    %             -1/4/pi * (Vel(XX)-Vel(YY))./vecnorm(XX-YY,2,2).^3;

    % Manually putting vel(x) = 0
    Ciekernel1 = 3/4/pi * (XX-YY).*dot(XX-YY,-Vel(YY),2)./vecnorm(XX-YY,2,2).^5 ...
                -1/4/pi * (-Vel(YY))./vecnorm(XX-YY,2,2).^3;

    Cie1 = sum(W.*dot(Psi_eXX,cross(Ciekernel1,nxgvals_iYY,2),2) ,1);
    
    sd = -1/(2*mu0) * ( 2 * (Aei1)...
                        + 2 * (Cie1));

end