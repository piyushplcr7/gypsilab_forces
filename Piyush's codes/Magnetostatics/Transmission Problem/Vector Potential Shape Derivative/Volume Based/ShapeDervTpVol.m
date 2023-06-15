% Shape derivative from volume based formulation
% Bn and Ht have to be evaluations at the X values for Gamma

function sd = ShapeDervTpVol(Gamma,Bn,Ht,mu_e,mu_i,vel)
    jump_mu_inv = 1/mu_e - 1/mu_i;
    jump_mu = mu_e - mu_i;

    [X,W] = Gamma.qud;
    normals = Gamma.qudNrm;
    vels = vel(X);

    sd = 0.5 * sum(W.* ((Bn).^2*jump_mu_inv - (Ht).^2*jump_mu).* dot(vels,normals,2), 1);
end