%
% Input Bn is the normal component of tilde B 
% Input Hot is the tangential tilde H field on the outside
% TdA is the dirichlet trace of the vector potential
% B0 is the constant field

function force = ForceMstTPConstFiend(Gamma,Bn,Hot,mu0,mu,B0,Vel,DVel,TdA)
    jump_mu_inv = 1/mu0 - 1/mu;
    jump_mu = mu0 - mu;

    bndmesh = Gamma.msh;

    [X,W] = Gamma.qud;
    normals = Gamma.qudNrm;

    Vels = Vel(X);
    DVel1XX = DVel{1}(X);
    DVel2XX = DVel{2}(X);
    DVel3XX = DVel{3}(X);

    % Constant field and its tangential component
    B_const_vals = B0(X);
    B_const_tau = cross(normals, cross(B_const_vals,normals,2) ,2);

    Hit = Hot + jump_mu_inv * B_const_tau;

    % Integrand for bilinear form b
    vectorintegrand = - 0.5 * jump_mu_inv * Bn.^2 .* normals...
                      + jump_mu_inv * Bn .* B_const_tau...
                      + 0.5 * normals .* (mu0 * vecnorm(Hot,2,2).^2 ...
                                        - mu * vecnorm(Hit,2,2).^2);
    vectorintegrand = dot(Vels,vectorintegrand,2);

    % Linear form
    NED = fem(bndmesh,'NED'); 
    g = reconstruct(TdA,Gamma,NED);
    nxA = cross(normals,g,2);
    integrandl = dot( [dot(DVel1XX,nxA,2)...
                       dot(DVel2XX,nxA,2)...
                       dot(DVel3XX,nxA,2)]...
                       ,B_const_vals,2);

    % Intuition term
    integrandintuition = ...
    -0.5 * jump_mu_inv * vecnorm(B_const_vals,2,2).^2.*dot(Vels,normals,2);
    integrandintuition = integrandintuition * 0;

%     contri1 = sum(W.* (vectorintegrand ) ,1)
%     contri2 = sum(W.* (- integrandl ) ,1)
%     contri3 = sum(W.* ( integrandintuition) ,1)

    force = sum(W.* (vectorintegrand - integrandl + integrandintuition) ,1);

end