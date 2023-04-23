function sd = shapeDer_ScalPot_Vol_TP(bndmesh,mu,mu0,psi,g,J,omega_src,Vel,DVel)
    % BEM spaces
    P1 = fem(bndmesh,'P1');
    P0 = fem(bndmesh,'P0');

    Gamma = dom(bndmesh,3);
    [X,Wts] = Gamma.qud;
    normals = Gamma.qudNrm;

    gvals = reconstruct(g,Gamma,P1);
    psivals = reconstruct(psi,Gamma,P0);
    gradTg = reconstruct(g,Gamma,P1.grad);

    HJ = compute_vecpot_curl(J,omega_src,X);
    HJn = dot(normals,HJ,2);
    HJn_coeffs = proj(HJn,Gamma,P0);

    DHJ = compute_vecpot_D_curl(J,omega_src,X);

    % Evaluating the velocity field at quadrature points
    Vels = Vel(X);
    % Evaluating the Jacobian of the velocity field (row-wise) at qud pts
    DVel1 = DVel{1}(X);
    DVel2 = DVel{2}(X);
    DVel3 = DVel{3}(X);
    divVel = DVel1(:,1) + DVel2(:,2) + DVel3(:,3);

    DHJV = [dot(DHJ{1},Vels,2) dot(DHJ{2},Vels,2) dot(DHJ{3},Vels,2)];
    DVHJ = [dot(DVel1,HJ,2) dot(DVel2,HJ,2) dot(DVel3,HJ,2)];

    jumpMu = mu0-mu;

    psi_in = mu0/mu*psi+jumpMu/mu*HJn_coeffs;
    psi_invals = reconstruct(psi_in,Gamma,P0);

    % t1
    % Tn_out
    Tno = 0.5*(psivals.^2-vecnorm(gradTg,2,2).^2) .* normals + psivals .*gradTg;
    % Tn_in
    Tni = 0.5*(psi_invals.^2-vecnorm(gradTg,2,2).^2) .* normals + psi_invals .*gradTg;
    jumpTn = mu0 * Tno - mu * Tni;
    t1 = sum(Wts.*dot(Vels,jumpTn,2),1);

    % t2
    t2 = -jumpMu * sum(Wts.*gvals.*dot(normals,(DHJV - DVHJ + HJ.*divVel),2),1);

    % t3
    Vn = dot(Vels,normals,2);
    t3 = -0.5 * jumpMu * sum(Wts.* dot(HJ,HJ,2) .* Vn,1);

    sd = t1+t2+t3;
end