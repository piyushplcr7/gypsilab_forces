% Maxwell stress tensor based evaluation
% Bn and Ht have to be evaluations at the X values for Gamma

function torque = TorqueMstTP(Gamma,Bn,Ht,mu_e,mu_i,Xcg)
    jump_mu_inv = 1/mu_e - 1/mu_i;
    jump_mu = mu_e - mu_i;

    [X,W] = Gamma.qud;
    normals = Gamma.qudNrm;
    rvec = X-Xcg;
    torque = 0.5 * sum(W.* ((Bn).^2*jump_mu_inv - (Ht).^2*jump_mu).* cross(rvec,normals,2), 1);
end