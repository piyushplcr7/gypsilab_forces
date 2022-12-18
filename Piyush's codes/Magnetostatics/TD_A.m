function trace = TD_A(J,omega_src,Gamma)
    % Need to compute the vector potential from the source Jsrc in
    % Omega_src and compute its dirichlet trace on Omega_target

    % Require Gamma_target to be a surface in 3 dim
    assert(size(Gamma.elt,2)==3);

    % Green's kernels
    Gxy = @(X,Y)1/4/pi*femGreenKernel(X,Y,'[1/r]',0); 

    % Obtaining quadrature nodes and weights
    [X_src,W_src,elt2qud_src] = omega_src.qud;
    [X,W,elt2qud] = Gamma.qud;
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

    % Taking the trace of the vector potential on Gamma
    normals = Gamma.msh.nrm;
    normals_X = repelem(normals,size(elt2qud,2),1);
    vecpot_X = [integral1 integral2 integral3];
    TD_vecpot_X = vecpot_X - dot(normals_X,vecpot_X,2).*normals_X;
    
    integral1 = TD_vecpot_X(:,1);
    integral2 = TD_vecpot_X(:,2);
    integral3 = TD_vecpot_X(:,3);

    % Projecting the traces to the space H^{-0.5}(curl_Gamma,Gamma)
    NED = fem(Gamma.msh,'NED');
    M = integral(Gamma,NED,NED);
    uqmat = NED.uqm(Gamma);
    % RHS is obtained by integrating TD_A.basis
    rhs = sum(W.*(uqmat{1}.*integral1+uqmat{2}.*integral2+uqmat{3}.*integral3),1)';
    trace = M\rhs;
    
end