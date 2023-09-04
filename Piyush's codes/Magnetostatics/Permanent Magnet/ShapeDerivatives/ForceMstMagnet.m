% Maxwell stress tensor based evaluation
% Bn and Ht have to be evaluations at the X values for Gamma. Ht passed as
% a vector

function force = ForceMstMagnet(Gamma,Bn,Ht,mu0,mu,M)
    jump_mu_inv = 1/mu0 - 1/mu;

    [X,W] = Gamma.qud;
    % Tangential component of magnetization
    Mvals = M(X);
    normals = Gamma.qudNrm;
    Mt = Mvals - dot(Mvals,normals,2).*normals;

    integrand = 0.5 * jump_mu_inv * Bn.^2 .* normals...
                - 0.5 * mu0 * vecnorm(Ht,2,2).^2 .* normals ...
                + 0.5 * mu * vecnorm(Ht+Mt,2,2).^2 .* normals ...
                - Mt.*Bn;

    force = sum(W.* integrand, 1);
end