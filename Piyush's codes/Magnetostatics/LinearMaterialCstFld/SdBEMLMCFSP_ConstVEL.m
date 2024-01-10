function sd = SdBEMLMCFSP_ConstVEL(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel,DVel,mu0,mu,H0)
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    jumpMu = mu0 - mu;

    % BEM Spaces
    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');
    P1_i = fem(bndmesh_i,'P1');
%     P1_e = fem(bndmesh_e,'P1');

    psi_i_vals = reconstruct(psi_i,Gamma_i,P0_i);
    g_i_vals = reconstruct(g_i,Gamma_i,P1_i);
    psi_e_vals = reconstruct(psi_e,Gamma_e,P0_e);

    normals_i = Gamma_i.qudNrm;

    %% Non SS Computation
    % partial derivatives Cross matrices
%     dsVei = single_layer_cross(Gamma_i,Gamma_e,P0_i,P0_e);
    dsVeiKernel = @(X,Y) 1/4/pi.* dot(Y-X,Vel(X),2) ./vecnorm(Y-X,2,2).^3;
    dsVei = integral(Gamma_i,Gamma_e,P0_i,dsVeiKernel,P0_e);

    % dsVei Explicit Computation
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

    % grady grady G(x,y) vel(y)
    dsKie1kernel = cell(3,1);
    dsKie1kernel{1} = @(X,Y) 3/4/pi * (X(:,1)-Y(1)) .* ((X-Y) * Vel(Y)')./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,Y,1) ./vecnorm(X-Y,2,2).^3;
    dsKie1kernel{2} = @(X,Y) 3/4/pi * (X(:,2)-Y(2)) .* ((X-Y) * Vel(Y)')./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,Y,2) ./vecnorm(X-Y,2,2).^3;
    dsKie1kernel{3} = @(X,Y) 3/4/pi * (X(:,3)-Y(3)) .* ((X-Y) * Vel(Y)')./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,Y,3) ./vecnorm(X-Y,2,2).^3;
    
    dsKie = integral(Gamma_e,Gamma_i,P0_e,dsKie1kernel,ntimes(P1_i));

    % dsKie explicit computation
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


    sd = -mu0/2*( 2 * psi_i' * dsVei * psi_e + 2 * psi_e' * dsKie * g_i);

%     pool.delete();
end