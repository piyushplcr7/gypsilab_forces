% Full shape Derivative

function sd = shapDervTranPrbScalPotBIE_CSTFLD(bndmesh,Tdu,Tnu,J,omega_src,Vel,DVel,mu0,mu)
    % BEM Spaces
    P0 = fem(bndmesh,'P0');
    P1 = fem(bndmesh,'P1');

    Gamma = dom(bndmesh,3);

    [Xgypsi,Wgypsi] = Gamma.qud;
    normals = Gamma.qudNrm;

    % Reconstructing the trace of u- > g
    g = reconstruct(Tdu,Gamma,P1);

    HJ = compute_vecpot_curl(J,omega_src,Xgypsi);
    DHJ = compute_vecpot_D_curl(J,omega_src,Xgypsi);

    HJn = dot(normals,HJ,2);
    HJn_coeffs = proj(HJn,Gamma,P0);

    % Evaluating the velocity field at quadrature points
    Vels = Vel(Xgypsi);
    % Evaluating the Jacobian of the velocity field (row-wise) at qud pts
    DVel1 = DVel{1}(Xgypsi);
    DVel2 = DVel{2}(Xgypsi);
    DVel3 = DVel{3}(Xgypsi);
    divVel = DVel1(:,1) + DVel2(:,2) + DVel3(:,3);

    jumpMu = mu0-mu;
    %% Full shape derivative

    %% Non SS computations

    % partial derivative of bk(g,psi)
    Kmat = double_layer_laplace(Gamma,P0,P1);

    %% SS Computations

    euler = parcluster('local');
    euler.NumWorkers = 5;
    saveProfile(euler);


    % Projecting the complex integrand to P0 space
    DHJV = [dot(DHJ{1},Vels,2) dot(DHJ{2},Vels,2) dot(DHJ{3},Vels,2)];
    DVHJ = [dot(DVel1,HJ,2) dot(DVel2,HJ,2) dot(DVel3,HJ,2)];
    compl_integrand = dot(normals,(DHJV - DVHJ + HJ.*divVel),2);
    compl_integrand_coeffs = proj(compl_integrand,Gamma,P0);
    slmat = single_layer(Gamma,P0,P0);

    % partial derivative of l1(psi)
    dl1_ds = Tnu' * (slmat * compl_integrand_coeffs);

    % partial derivative of l2(g)
    dl2_ds = sum(Wgypsi.*compl_integrand.*g,1);

    % partial derivative of l3(g)
    dl3_ds = compl_integrand_coeffs' * Kmat * Tdu;
    
    % Remaining terms
    Vn = dot(Vels,normals,2);
    remterm1 = -0.5 * jumpMu * sum(Wgypsi.* dot(HJ,HJ,2) .* Vn,1);
    remterms = -jumpMu^2/2/mu *HJn_coeffs' * (2* slmat * compl_integrand_coeffs);

    sd = -mu0 * ( jumpMu/mu * dl1_ds ...
                  + jumpMu/2/mu0 * dl2_ds ...
                  -jumpMu/mu0 * dl3_ds)...
                  +remterm1 + remterms;

    % fileID = fopen('contributions.txt','a');
    % fprintf(fileID,'%.16f  %.16f  %.16f  %.16f  %.16f\n', full(dl1_ds), full(dl2_ds), full(dl3_ds), full(remterm1), full(remterms));

    % disp('printing shape derivative contributions');
    % contributions = [dl1_ds dl2_ds dl3_ds remterm1 remterms];
    % contributions = full(contributions)

end