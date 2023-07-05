function sd = PermanentMagnetShapeDerivativeBIESP_BACKUP(Gamma,Tdu,Tnu,J,omega_src,Vel,DVel,mu0,M)
    % BEM Spaces
    bndmesh = Gamma.msh;
    P0 = fem(bndmesh,'P0');
    P1 = fem(bndmesh,'P1');
    gradP1 = P1.grad;
    nxgradP1 = P1.nxgrad;

    Gamma = dom(bndmesh,3);

    [X,W] = Gamma.qud;
    normals = Gamma.qudNrm;

    % Reconstructing the Neumann trace psi
    psi = reconstruct(Tnu,Gamma,P0);
    % Reconstructing the surface gradient of u
    gradTu = reconstruct(Tdu,Gamma,gradP1);
    % Reconstructing the trace of u- > g
    g = reconstruct(Tdu,Gamma,P1);

    HJ = compute_vecpot_curl(J,omega_src,X);
    DHJ = compute_vecpot_D_curl(J,omega_src,X);

    % Evaluating the velocity field at quadrature points
    Vels = Vel(X);
    % Evaluating the Jacobian of the velocity field (row-wise) at qud pts
    DVel1 = DVel{1}(X);
    DVel2 = DVel{2}(X);
    DVel3 = DVel{3}(X);
    divVel = DVel1(:,1) + DVel2(:,2) + DVel3(:,3);

    %% bilinear form
    Nelt = bndmesh.nelt;

    [ii,jj] = meshgrid(1:Nelt,1:Nelt);

    % Kernel gradxG.vel(x) + gradyG.vel(y), z:= y-x
    kernelold = @(x,y,z) dot(z,Vel(x) - Vel(y), 2)./(vecnorm(z,2,2).^3)/ (4*pi);

    kerneloldmat_P0_P0 = panel_assembly(bndmesh,kernelold,P0,P0,ii(:),jj(:));
    kerneloldmat_nxgradP1_nxgradP1 = panel_assembly(bndmesh,kernelold,nxgradP1,nxgradP1,ii(:),jj(:));

    % Partial derivative of bv(psi,psi)
    dbv_ds = Tnu' * kerneloldmat_P0_P0 * Tnu;

    % partial derivative of bk(g,psi)
    divVelg = divVel.*g;
    divVelg_coeffs = proj(divVelg,Gamma,P1);
    Kmat = double_layer_laplace(Gamma,P0,P1);

    kernelintegrable = @(x,y,z) 3/(4*pi)* dot(z,Vel(y) - Vel(x),2) .*z ./vecnorm(z,2,2).^5;
    kernelintegrablemat = panel_assembly(bndmesh,kernelintegrable,ntimes(P1),P0,ii(:),jj(:));

    % Combination kernel that cancels singularity
    combkernel = @(x,y,z) 1/(4*pi) * ( -[ dot(DVel{1}(y),z,2) dot(DVel{2}(y),z,2) dot(DVel{3}(y),z,2) ] + Vel(y) - Vel(x) )./vecnorm(z,2,2).^3;
    combkernelmat = panel_assembly(bndmesh,combkernel,ntimes(P1),P0,ii(:),jj(:));

    dbk_ds = Tnu' * Kmat * divVelg_coeffs + Tnu' * (kernelintegrablemat -combkernelmat) * Tdu;

    % Partial derivative of bw(g,g)
    KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
    SL_Dvelnxgrad_nxgrad = panel_assembly_shape_derivative(bndmesh,KV,nxgradP1,nxgradP1,ii(:),jj(:),Vel,DVel);
    dbw_ds = Tdu' * ( kerneloldmat_nxgradP1_nxgradP1 + 2 * SL_Dvelnxgrad_nxgrad) * Tdu;

    %% Linear form
    Mvals = M(X);
    Mdotn = dot(Mvals,normals,2);

    lambda_coeffs = proj(Mdotn,Gamma,P0);

    l1 = Tnu' * kerneloldmat_P0_P0 * lambda_coeffs;

    l2 = lambda_coeffs' * Kmat * divVelg_coeffs + lambda_coeffs' * (kernelintegrablemat -combkernelmat) * Tdu;
    l2 = -l2;

    %% Remaining terms

    r1 = -mu0/2 * (lambda_coeffs' * kerneloldmat_P0_P0 * lambda_coeffs);
    r2 = mu0 * sum(W.* Mdotn .* dot(Vels,HJ,2),1);


    sd = mu0/2 * (-2*dbv_ds + 4 * dbk_ds +2 * dbw_ds)+ mu0  * (l1 + l2) +r1 + r2;
    

end