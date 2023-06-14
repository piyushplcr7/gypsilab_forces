% Compute the shape derivative obtained from the volume based scalar
% potential formulation of the superconductor

% Super conductor Shape derivative_ Scalar Potential
function sd = ScSd_SP_BEM_BACKUP(bndmesh,Tdu,Tnu,J,omega_src,Vel,DVel)
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

    %% Evaluating the three terms in the shape derivative

    % Evaluating the velocity field at quadrature points
    Vels = Vel(X);
    % Evaluating the Jacobian of the velocity field (row-wise) at qud pts
    DVel1 = DVel{1}(X);
    DVel2 = DVel{2}(X);
    DVel3 = DVel{3}(X);
    divVel = DVel1(:,1) + DVel2(:,2) + DVel3(:,3);

%     % Stress tensor times normal: T(u) n
%     Tun = 0.5 * (psi.^2 - dot(gradTu,gradTu,2)) .* normals + psi.* gradTu;
%     t1 = sum(W.* dot(Vels, Tun,2) ,1);
% 
%     % 3rd term
%     Vn = dot(Vels,normals,2);
%     t3 = -0.5 * sum(W.* dot(HJ,HJ,2) .* Vn,1);
% 
%     % 2nd term
%     DHJV = [dot(DHJ{1},Vels,2) dot(DHJ{2},Vels,2) dot(DHJ{3},Vels,2)];
%     DVHJ = [dot(DVel1,HJ,2) dot(DVel2,HJ,2) dot(DVel3,HJ,2)];
%     t2 = -sum(W.*g.*dot(normals,(DHJV - DVHJ + HJ.*divVel),2),1);
% 
%     sd = t1+t2+t3;

    %% Full shape derivative
    Nelt = bndmesh.nelt;

    [ii,jj] = meshgrid(1:Nelt,1:Nelt);

    % t1
    % Kernel for t1, z:= y-x
    kernelt1 = @(x,y,z) sum(z.*(Vel(x) - Vel(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    t1mat = panel_assembly(bndmesh,kernelt1,nxgradP1,nxgradP1,ii(:),jj(:));
    t1 = 0.5 * Tdu' * t1mat * Tdu;

    % t2
    KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
    t2mat = panel_assembly_shape_derivative(bndmesh,KV,nxgradP1,nxgradP1,ii(:),jj(:),Vel,DVel);
    t2 = Tdu' * t2mat * Tdu;

    % t3
    % Evaluating the velocity field at quadrature points
    Vels = Vel(X);
    % Evaluating the Jacobian of the velocity field (row-wise) at qud pts
    DVel1 = DVel{1}(X);
    DVel2 = DVel{2}(X);
    DVel3 = DVel{3}(X);
    divVel = DVel1(:,1) + DVel2(:,2) + DVel3(:,3);
    
    DHJV = [dot(DHJ{1},Vels,2) dot(DHJ{2},Vels,2) dot(DHJ{3},Vels,2)];
    DVHJ = [dot(DVel1,HJ,2) dot(DVel2,HJ,2) dot(DVel3,HJ,2)];
    t3 = -0.5*sum(W.*g.*dot(normals,(DHJV - DVHJ + HJ.*divVel),2),1);

    % t4 
    % Projecting the complex integrand to P0 space
    compl_integrand = dot(normals,(DHJV - DVHJ + HJ.*divVel),2);
    compl_integrand_coeffs = proj(compl_integrand,Gamma,P0);
    Kmat = double_layer_laplace(Gamma,P0,P1);
    t4 = -compl_integrand_coeffs' * Kmat * Tdu;

    % t8
    HJn = dot(normals,HJ,2);
    HJn_coeffs = proj(HJn,Gamma,P0);
    Vmat = single_layer(Gamma,P0,P0);
    t8 = -HJn_coeffs' * Vmat * compl_integrand_coeffs;

    % t7
    t7mat = panel_assembly(bndmesh,kernelt1,P0,P0,ii(:),jj(:));
    t7 = -0.5 * HJn_coeffs' * t7mat * HJn_coeffs;

    % t9
    Vn = dot(Vels,normals,2);
    t9 = -0.5 * sum(W.* dot(HJ,HJ,2) .* Vn,1);

    % t5 and t6
    % Double layer like term in t5
    divVelg = divVel.*g;
    divVelg_coeffs = proj(divVelg,Gamma,P1);
    t5dl = -HJn_coeffs' *Kmat * divVelg_coeffs;
    % integrable kernel in t6
    kernelt6 = @(x,y,z) 3/(4*pi)* dot(z,Vel(y) - Vel(x),2) .*z ./vecnorm(z,2,2).^5;
    kernelt6mat = panel_assembly(bndmesh,kernelt6,ntimes(P1),P0,ii(:),jj(:));
    t6 = -HJn_coeffs' * kernelt6mat * Tdu;
    % Combination kernel of t5 and t6 that cancels singularity
    combkernel = @(x,y,z) 1/(4*pi) * ( -[ dot(DVel{1}(y),z,2) dot(DVel{2}(y),z,2) dot(DVel{3}(y),z,2) ] + Vel(y) - Vel(x) )./vecnorm(z,2,2).^3;
    combkernelmat = panel_assembly(bndmesh,combkernel,ntimes(P1),P0,ii(:),jj(:));
    t56 = HJn_coeffs' * combkernelmat * Tdu;

    sd = t1+t2+t3+t4+t5dl+t56+t6+t7+t8+t9;

end