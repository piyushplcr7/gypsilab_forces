% Computes the quantity
%   \int_{\Omega_src} G(x,y) J(y) dy
%   the points x are specified in the input X

function A = compute_vecpot(J,omega_src,X)
    % Green's kernels
    Gxy = @(X,Y)1/4/pi*femGreenKernel(X,Y,'[1/r]',0); 
    [X_src,W_src,~] = omega_src.qud;
    N_qud_pts = size(X,1);
    N_qud_pts_src = size(X_src,1);

    % Evaluating source current at source quadrature nodes
    J_X_src = J(X_src(:,1),X_src(:,2),X_src(:,3));

    % quadrature entries to be summed without greens function
    integrand_without_green_fn = W_src.*J_X_src;

    % Contains [X_src; X_src; ... ]
    X_src_extended = repmat(X_src,N_qud_pts,1);
    % Contains [X(1,:); X(1,:)......; X(2,:); X(2,:); ....]
    X_extended = repelem(X,N_qud_pts_src,1);

    % Evaluating the green's function
    green_fn = Gxy(X_src_extended,X_extended);
    green_fn = reshape(green_fn,[N_qud_pts_src, N_qud_pts]);

    % Vector potential evaluated at X
    integral1 = sum(integrand_without_green_fn(:,1).*green_fn,1)';
    integral2 = sum(integrand_without_green_fn(:,2).*green_fn,1)';
    integral3 = sum(integrand_without_green_fn(:,3).*green_fn,1)';

    A = [integral1 integral2 integral3];

end