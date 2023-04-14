% Vectorial Single Layer Potential using coefficients of a density in
% the discrete space RWG,evaluated at the points in X (NX3)

% Computes the quantity
%   \int_{\Gamma} G(x,y) psi(y) dy
%   the points x are specified in the input X

function A = SLP_Vectorial(Gamma,psi_coeffs,X)
    % Green's kernel
    Gxy = @(X,Y)1/4/pi*femGreenKernel(X,Y,'[1/r]',0); 
    [X_Gamma,W_Gamma,~] = Gamma.qud;
    N_eval_pts = size(X,1);
    N_qud_pts_Gamma = size(X_Gamma,1);

    % Evaluating source density at source quadrature nodes
    % BEM space in which the source density lies
    RWG = fem(Gamma.msh,'RWG');
    psi = reconstruct(psi_coeffs,Gamma,RWG);

    % quadrature entries to be summed without greens function
    integrand_without_green_fn = W_Gamma.*psi;

    % Contains [X_Gamma; X_Gamma; ... ]
    X_Gamma_extended = repmat(X_Gamma,N_eval_pts,1);
    % Contains [X(1,:); X(1,:)......; X(2,:); X(2,:); ....]
    X_extended = repelem(X,N_qud_pts_Gamma,1);

    % Evaluating the green's function
    green_fn = Gxy(X_Gamma_extended,X_extended);
    green_fn = reshape(green_fn,[N_qud_pts_Gamma, N_eval_pts]);

    % Vector potential evaluated at X
    integral1 = sum(integrand_without_green_fn(:,1).*green_fn,1)';
    integral2 = sum(integrand_without_green_fn(:,2).*green_fn,1)';
    integral3 = sum(integrand_without_green_fn(:,3).*green_fn,1)';

    A = [integral1 integral2 integral3];

end