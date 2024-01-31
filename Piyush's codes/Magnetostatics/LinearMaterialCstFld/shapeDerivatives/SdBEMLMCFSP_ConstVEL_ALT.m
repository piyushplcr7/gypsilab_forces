function sd = SdBEMLMCFSP_ConstVEL_ALT(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel,DVel,mu0,mu,H0)
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    % jumpMu = mu0 - mu;

    % BEM Spaces
    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');
    P1_i = fem(bndmesh_i,'P1');
    P1_e = fem(bndmesh_e,'P1');

    psi_i_vals = reconstruct(psi_i,Gamma_i,P0_i);
    g_i_vals = reconstruct(g_i,Gamma_i,P1_i);
    grad_g_i_vals = reconstruct(g_i,Gamma_i,P1_i.grad);
    psi_e_vals = reconstruct(psi_e,Gamma_e,P0_e);

    % Dirichlet trace at the outer boundary
    [X_e,~] = Gamma_e.qud;
    g_e_vals = X_e * H0';
    g_e_coeffs = proj(g_e_vals,Gamma_e,P1_e);

    normals_i = Gamma_i.qudNrm;

    %% Non SS Computation
    % partial derivatives Cross matrices

    % gradx G(x,y) . vel(x)
    dsVeiKernel = @(X,Y) 1/4/pi.* dot(Y-X,Vel(X),2) ./vecnorm(Y-X,2,2).^3;
    dsVei = integral(Gamma_i,Gamma_e,P0_i,dsVeiKernel,P0_e);

    % % % dsVei Explicit Computation
    % [X_i,W_i] = Gamma_i.qud;
    % [Y_e,W_e] = Gamma_e.qud;
    % NX = size(X_i,1);
    % NY = size(Y_e,1);
    % 
    % XX = repelem(X_i,NY,1); psi_i_XX = repelem(psi_i_vals,NY,1);
    % YY = repmat(Y_e,NX,1); psi_e_YY = repmat(psi_e_vals,NX,1);
    % 
    % W = repelem(W_i,NY,1).*repmat(W_e,NX,1);
    % 
    % kernel = 1/4/pi * dot(YY-XX,Vel(XX),2)./vecnorm(XX-YY,2,2).^3;
    % 
    % testyo = sum(W.*kernel.*psi_e_YY.*psi_i_XX ,1);
    % 
    % % dsl1 explicit computation
    % H0dotYY = YY * H0';
    % normals_e = Gamma_e.qudNrm;
    % normals_e_YY = repmat(normals_e,NX,1);
    % kdsl1 = 1/4/pi * Vel(XX)./vecnorm(XX-YY,2,2).^3 - 3/4/pi * (XX-YY).*dot(XX-YY,Vel(XX),2)./vecnorm(XX-YY,2,2).^5;
    % dsl1explicit = sum(W.*dot(kdsl1,normals_e_YY,2).*H0dotYY.*psi_i_XX,1);
    % 
    % % dsl2 explicit computation
    % H0_YY = repmat(H0,size(YY,1),1);
    % nxH0YY = cross(normals_e_YY,H0_YY,2);
    % gradg_i_vals_XX = repelem(grad_g_i_vals,NY,1);
    % normals_i_XX = repelem(normals_i,NY,1);
    % kdsl2 = 1/4/pi * dot(XX-YY,Vel(XX),2)./vecnorm(XX-YY,2,2).^3;
    % dsl2explicit = sum(W.*kdsl2.*dot(nxH0YY,cross(normals_i_XX,gradg_i_vals_XX,2),2),1);

    % grady grady G(x,y) vel(y)
    dsKie1kernel = cell(3,1);
    dsKie1kernel{1} = @(X,Y) 3/4/pi * (X(:,1)-Y(1)) .* ((X-Y) * Vel(Y)')./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,Y,1) ./vecnorm(X-Y,2,2).^3;
    dsKie1kernel{2} = @(X,Y) 3/4/pi * (X(:,2)-Y(2)) .* ((X-Y) * Vel(Y)')./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,Y,2) ./vecnorm(X-Y,2,2).^3;
    dsKie1kernel{3} = @(X,Y) 3/4/pi * (X(:,3)-Y(3)) .* ((X-Y) * Vel(Y)')./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,Y,3) ./vecnorm(X-Y,2,2).^3;
    
    dsKie = integral(Gamma_e,Gamma_i,P0_e,dsKie1kernel,ntimes(P1_i));

    dsl1kernel = cell(3,1); % grad x grad x G vel(x)
    dsl1kernel{1} = @(X,Y) 3/4/pi * (X(:,1)-Y(1)) .* dot((X-Y), Vel(X),2)./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,X,1) ./vecnorm(X-Y,2,2).^3;
    dsl1kernel{2} = @(X,Y) 3/4/pi * (X(:,2)-Y(2)) .* dot((X-Y), Vel(X),2)./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,X,2) ./vecnorm(X-Y,2,2).^3;
    dsl1kernel{3} = @(X,Y) 3/4/pi * (X(:,3)-Y(3)) .* dot((X-Y), Vel(X),2)./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,X,3) ./vecnorm(X-Y,2,2).^3;
    dsl1 = -integral(Gamma_i,Gamma_e,P0_i,dsl1kernel,ntimes(P1_e));


    dsl2 = -integral(Gamma_i,Gamma_e,P1_i.nxgrad,dsVeiKernel,P1_e.nxgrad);

    % % % dsKie explicit computation
    % [Y,WY] = Gamma_i.qud;
    % [X,WX] = Gamma_e.qud;
    % 
    % NX = size(X,1); NY = size(Y,1);
    % 
    % XX = repelem(X,NY,1); YY = repmat(Y,NX,1);
    % W = repelem(WX,NY,1).*repmat(WY,NX,1);
    % 
    % g_i_vals_YY = repmat(g_i_vals,NX,1);
    % psi_e_XX = repelem(psi_e_vals,NY,1);
    % normals_yy = repmat(normals_i,NX,1);
    % 
    % kernel = 3/4/pi * (XX-YY).*dot(XX-YY,Vel(YY),2)./vecnorm(XX-YY,2,2).^5 - 1/4/pi * Vel(YY)./vecnorm(XX-YY,2,2).^3;
    % 
    % testa = sum(W.* dot(kernel,normals_yy,2) .* g_i_vals_YY.*psi_e_XX ,1);


    sd = -mu0/2*( 2 * psi_i' * dsVei * psi_e + 2 * psi_e' * dsKie * g_i)...
        + mu0 * (psi_i' * dsl1 * g_e_coeffs + g_i' * dsl2 * g_e_coeffs);

%     pool.delete();
end