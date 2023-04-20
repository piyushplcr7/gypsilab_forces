% Compute the shape derivative obtained from the volume based scalar
% potential formulation of the superconductor

% Super conductor Shape derivative_ Scalar Potential
function sd = ScSd_SP_BEM(bndmesh,Tdu,Tnu,J,omega_src,Vel,DVel)
    % BEM Spaces
    P0 = fem(bndmesh,'P0');
    P1 = fem(bndmesh,'P1');
    gradP1 = P1.grad;

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

    %% Evaluating the three terms in the shape derivative

    % Evaluating the velocity field at quadrature points
    Vels = Vel(X);
    % Evaluating the Jacobian of the velocity field (row-wise) at qud pts
    DVel1 = DVel{1}(X);
    DVel2 = DVel{2}(X);
    DVel3 = DVel{3}(X);
    divVel = DVel1(:,1) + DVel2(:,2) + DVel3(:,3);

    % Stress tensor times normal: T(u) n
    Tun = 0.5 * (psi.^2 - dot(gradTu,gradTu,2)) .* normals + psi.* gradTu;
    t1 = sum(W.* dot(Vels, Tun,2) ,1);

    % 3rd term
    Vn = dot(Vels,normals,2);
    t3 = -0.5 * sum(W.* dot(HJ,HJ,2) .* Vn,1);

    % 2nd term
    DHJV = [dot(DHJ{1},Vels,2) dot(DHJ{2},Vels,2) dot(DHJ{3},Vels,2)];
    DVHJ = [dot(DVel1,HJ,2) dot(DVel2,HJ,2) dot(DVel3,HJ,2)];
    t2 = -sum(W.*g.*dot(normals,(DHJV - DVHJ + HJ.*divVel),2),1);

    sd = t1+t2+t3;

end