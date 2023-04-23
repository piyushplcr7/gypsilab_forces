% Maxwell stress tensor based evaluation
% Bn and Ht have to be evaluations at the X values for Gamma

function force = ForceMstTP(Gamma,Bn,Ht,mu_e,mu_i)
    jump_mu_inv = 1/mu_e - 1/mu_i;
    jump_mu = mu_e - mu_i;

    [X,W] = Gamma.qud;
    normals = Gamma.qudNrm;
    force = 0.5 * sum(W.* ((Bn).^2*jump_mu_inv - (Ht).^2*jump_mu).* normals, 1);
end