% Compute the shape derivative obtained from the BEM based scalar
% potential formulation of the superconductor. Simplified for constant
% Fields which are used for computing forces

% Super conductor Shape derivative_ Scalar Potential
function sd = ScSdF_SP_BEM(bndmesh,Tdu,Tnu,J,omega_src,Vel)
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

    %% Evaluating the shape derivative for constant fields

    % Evaluating the velocity field at quadrature points
    Vels = Vel(X);
    DHJV = [dot(DHJ{1},Vels,2) dot(DHJ{2},Vels,2) dot(DHJ{3},Vels,2)];

    % t1
    t1 = -0.5 * sum(W.*g.*dot(normals,DHJV,2),1);

    % t2 Computing t2 by projecting the x - integrand to P0 first
    nT_DHJV = dot(normals,DHJV,2);
    nT_DHJV_coeffs = proj(nT_DHJV,Gamma,P0);

    K = double_layer_laplace(Gamma,P0,P1);
    t2 = -nT_DHJV_coeffs' * K * Tdu;

    % 3rd term
    Vn = dot(Vels,normals,2);
    t3 = -0.5 * sum(W.* dot(HJ,HJ,2) .* Vn,1);

    % 4th term
    HJn = dot(normals,HJ,2);
    HJn_coeffs = proj(HJn,Gamma,P0);
    V = single_layer(Gamma,P0,P0);
    t4 = -HJn_coeffs' * V * nT_DHJV_coeffs;

    
    sd = t1 + t2 + t3 + t4;

end