function sd = SdBEMLMCFSP_ALT(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel,DVel,mu0,mu,H0)
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    jumpMu = mu0 - mu;

    % BEM Spaces
    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');
    P1_i = fem(bndmesh_i,'P1');
    P1_e = fem(bndmesh_e,'P1');

    psi_i_vals = reconstruct(psi_i,Gamma_i,P0_i);
    g_i_vals = reconstruct(g_i,Gamma_i,P1_i);
    psi_e_vals = reconstruct(psi_e,Gamma_e,P0_e);
    normals_i = Gamma_i.qudNrm;
    normals_e = Gamma_e.qudNrm;
    grad_g_i_vals = reconstruct(g_i,Gamma_i,P1_i.grad);

    % Dirichlet trace at the outer boundary
    [X_e,~] = Gamma_e.qud;
    g_e_vals = X_e * H0';
    g_e_coeffs = proj(g_e_vals,Gamma_e,P1_e);

    %% Non SS Computations
    % Evaluating the Jacobian of the velocity field (row-wise) at gamma_i qud pts
    [X_i,W_i] = Gamma_i.qud;
    Veli = Vel(X_i);
    DVel1i = DVel{1}(X_i);
    DVel2i = DVel{2}(X_i);
    DVel3i = DVel{3}(X_i);
    divVeli = DVel1i(:,1) + DVel2i(:,2) + DVel3i(:,3);
    divVelgi = divVeli.*g_i_vals;
    divVelgi_coeffs = proj(divVelgi,Gamma_i,P1_i);

    Kmatii = double_layer_laplace(Gamma_i,P0_i,P1_i);
    Vii = single_layer(Gamma_i,P0_i,P0_i);

    % partial derivatives Cross matrices
    % gradx G(x,y) . vel(x)
    dsVeiKernel = @(X,Y) 1/4/pi.* dot(Y-X,Vel(X),2) ./vecnorm(Y-X,2,2).^3;
    dsVei = integral(Gamma_i,Gamma_e,P0_i,dsVeiKernel,P0_e);

    % grady grady G(x,y) vel(y)
    dsKie1kernel = cell(3,1);
    dsKie1kernel{1} = @(X,Y) 3/4/pi * (X(:,1)-Y(1)) .* ((X-Y) * Vel(Y)')./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,Y,1) ./vecnorm(X-Y,2,2).^3;
    dsKie1kernel{2} = @(X,Y) 3/4/pi * (X(:,2)-Y(2)) .* ((X-Y) * Vel(Y)')./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,Y,2) ./vecnorm(X-Y,2,2).^3;
    dsKie1kernel{3} = @(X,Y) 3/4/pi * (X(:,3)-Y(3)) .* ((X-Y) * Vel(Y)')./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,Y,3) ./vecnorm(X-Y,2,2).^3;
    dsKie1 = integral(Gamma_e,Gamma_i,P0_e,dsKie1kernel,ntimes(P1_i));

    % dsKie2 explicit computation
    [Y,WY] = Gamma_i.qud;
    [X,WX] = Gamma_e.qud;

    NX = size(X,1); NY = size(Y,1);

    XX = repelem(X,NY,1); YY = repmat(Y,NX,1);
    W = repelem(WX,NY,1).*repmat(WY,NX,1);

    g_i_vals_YY = repmat(g_i_vals,NX,1);
    psi_e_XX = repelem(psi_e_vals,NY,1);
    normals_yy = repmat(normals_i,NX,1);
    DVel1YY = DVel{1}(YY);
    DVel2YY = DVel{2}(YY);
    DVel3YY = DVel{3}(YY);
    divVelYY = DVel1YY(:,1)+DVel2YY(:,2)+DVel3YY(:,3);
    dsKie2kernel =1/4/pi./vecnorm(XX-YY).^3  .*(divVelYY.*(XX-YY) - [dot(DVel1YY,XX-YY,2) dot(DVel2YY,XX-YY,2) dot(DVel3YY,XX-YY,2)]);
    dsKie2 = sum(W.* dot(dsKie2kernel,normals_yy,2) .* g_i_vals_YY.*psi_e_XX ,1);
    
    % dsl3 explicit computation
    kerneldsl3 = 1/4/pi./vecnorm(XX-YY,2,2);
    % grad_g_i_vals_YY = repmat(grad_g_i_vals,NX,1);
    normals_e_XX = repelem(normals_e,NY,1);
    H0XX = repelem(H0,size(normals_e_XX,1),1);
    nxH0XX = cross(normals_e_XX,H0XX,2);
    nxgradgiYY = repmat(cross(normals_i,grad_g_i_vals,2),NX,1);
    DVelnxgradgiYY = [dot(DVel1YY,nxgradgiYY,2) dot(DVel2YY,nxgradgiYY,2) dot(DVel3YY,nxgradgiYY,2)];
    dsl3 = -sum(W.*kerneldsl3.*dot(DVelnxgradgiYY,nxH0XX,2),1);

    % Vee = single_layer(Gamma_e,P0_e,P0_e);

    % H0.n coefficients
    % H0extended = repmat(H0,size(normals_i,1),1);
    % H0dotn_vals = dot(H0extended,normals_i,2);
    % H0dotn_coeffs = proj(H0dotn_vals,Gamma_i,P0_i);
    % 
    % % Complex integrand
    % DVelH0 = [dot(DVel1i,H0extended,2) dot(DVel2i,H0extended,2) dot(DVel3i,H0extended,2)];
    % compl_integrand = dot(normals_i,H0extended.*divVeli - DVelH0,2);
    % compl_integrand_coeffs = proj(compl_integrand,Gamma_i,P0_i);

    %% SS Computations
    Nelt_i = bndmesh_i.nelt;

    [ii,jj] = meshgrid(1:Nelt_i,1:Nelt_i);

    % Kernel gradxG.vel(x) + gradyG.vel(y), z:= y-x
    kernelold = @(x,y,z) dot(z,Vel(x) - Vel(y), 2)./(vecnorm(z,2,2).^3)/ (4*pi);

    kernelintegrable = @(x,y,z) 3/(4*pi)* dot(z,Vel(y) - Vel(x),2) .*z ./vecnorm(z,2,2).^5;

    combkernel = @(x,y,z) 1/(4*pi) * ( -[ dot(DVel{1}(y),z,2) dot(DVel{2}(y),z,2) dot(DVel{3}(y),z,2) ] + Vel(y) - Vel(x) )./vecnorm(z,2,2).^3;
    
    KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;

    euler = parcluster('local');
    euler.NumWorkers = 5;
    saveProfile(euler);

    pool = euler.parpool(5);

    spmd
        if spmdIndex==1
            kerneloldmat_P0_P0_ii = panel_assembly(bndmesh_i,kernelold,P0_i,P0_i,ii(:),jj(:));
            % Partial derivative of bv(psi,psi)
            dbv_dsii = psi_i' * kerneloldmat_P0_P0_ii * psi_i;

        elseif spmdIndex==2
            kerneloldmat_nxgradP1_nxgradP1ii = panel_assembly(bndmesh_i,kernelold,P1_i.nxgrad,P1_i.nxgrad,ii(:),jj(:));

        elseif spmdIndex==3
            kernelintegrablematii = panel_assembly(bndmesh_i,kernelintegrable,ntimes(P1_i),P0_i,ii(:),jj(:));

        elseif spmdIndex==4
            % Combination kernel that cancels singularity
            combkernelmatii = panel_assembly(bndmesh_i,combkernel,ntimes(P1_i),P0_i,ii(:),jj(:));

        elseif spmdIndex==5
             % Partial derivative of bw(g,g)
            SL_Dvelnxgrad_nxgradii = panel_assembly_shape_derivative(bndmesh_i,KV,P1_i.nxgrad,P1_i.nxgrad,ii(:),jj(:),Vel,DVel);
        end
    end

    dbk_dsii = psi_i' * Kmatii * divVelgi_coeffs + psi_i' * (kernelintegrablematii{3} -combkernelmatii{4}) * g_i;
    dbw_dsii = g_i' * ( kerneloldmat_nxgradP1_nxgradP1ii{2} + 2 * SL_Dvelnxgrad_nxgradii{5}) * g_i;

    % Linear form
    dsl1kernel = cell(3,1); % grad x grad x G vel(x)
    dsl1kernel{1} = @(X,Y) 3/4/pi * (X(:,1)-Y(1)) .* dot((X-Y), Vel(X),2)./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,X,1) ./vecnorm(X-Y,2,2).^3;
    dsl1kernel{2} = @(X,Y) 3/4/pi * (X(:,2)-Y(2)) .* dot((X-Y), Vel(X),2)./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,X,2) ./vecnorm(X-Y,2,2).^3;
    dsl1kernel{3} = @(X,Y) 3/4/pi * (X(:,3)-Y(3)) .* dot((X-Y), Vel(X),2)./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,X,3) ./vecnorm(X-Y,2,2).^3;
    dsl1 = -integral(Gamma_i,Gamma_e,P0_i,dsl1kernel,ntimes(P1_e));


    dsl2 = -integral(Gamma_i,Gamma_e,P1_i.nxgrad,dsVeiKernel,P1_e.nxgrad);


    sd = -mu0/2*( (1+mu0/mu) * dbv_dsii{1} + 4 * dbk_dsii -(1+mu/mu0) * dbw_dsii...
                  + 2 * psi_i' * dsVei * psi_e + 2 * psi_e' * dsKie1 * g_i + 2 * dsKie2)...
         +mu0 * (psi_i' * dsl1 * g_e_coeffs + g_i' * dsl2 * g_e_coeffs + dsl3);

    pool.delete();
end