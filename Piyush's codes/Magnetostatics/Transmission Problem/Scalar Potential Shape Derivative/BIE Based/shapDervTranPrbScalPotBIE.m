% Full shape Derivative

function sd = shapDervTranPrbScalPotBIE(bndmesh,Tdu,Tnu,J,omega_src,Vel,DVel,mu0,mu)
    % BEM Spaces
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

    jumpMu = mu0-mu;
    %% Full shape derivative
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

    % Projecting the complex integrand to P0 space
    DHJV = [dot(DHJ{1},Vels,2) dot(DHJ{2},Vels,2) dot(DHJ{3},Vels,2)];
    DVHJ = [dot(DVel1,HJ,2) dot(DVel2,HJ,2) dot(DVel3,HJ,2)];
    compl_integrand = dot(normals,(DHJV - DVHJ + HJ.*divVel),2);
    compl_integrand_coeffs = proj(compl_integrand,Gamma,P0);
    HJn = dot(normals,HJ,2);
    HJn_coeffs = proj(HJn,Gamma,P0);
    slmat = single_layer(Gamma,P0,P0);

    % partial derivative of l1(psi)
    dl1_ds = Tnu' * (kerneloldmat_P0_P0 * HJn_coeffs + slmat * compl_integrand_coeffs);

    % partial derivative of l2(g)
    dl2_ds = sum(W.*compl_integrand.*g,1);

    % partial derivative of l3(g)
    dl3_ds = HJn_coeffs' * Kmat * divVelg_coeffs + HJn_coeffs' * (kernelintegrablemat -combkernelmat) * Tdu...
            + compl_integrand_coeffs' * Kmat * Tdu;
    
    % Remaining terms
    Vn = dot(Vels,normals,2);
    remterm1 = -0.5 * jumpMu * sum(W.* dot(HJ,HJ,2) .* Vn,1);
    remterms = -jumpMu^2/2/mu *HJn_coeffs' * (kerneloldmat_P0_P0 * HJn_coeffs + 2* slmat * compl_integrand_coeffs);

    sd = mu0/2 * ( - (1+mu0/mu) * dbv_ds ...
                   + 4 * dbk_ds...
                   + (1+mu/mu0) * dbw_ds)...
         -mu0 * ( jumpMu/mu * dl1_ds ...
                  + jumpMu/2/mu0 * dl2_ds ...
                  -jumpMu/mu0 * dl3_ds)...
                  +remterm1 + remterms;



end