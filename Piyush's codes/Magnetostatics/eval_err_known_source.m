function [l2err,l2norm_sol,l2norm] = eval_err_known_source(N)

    L = [1 1 1];
    mesh = mshCube(N,L);
    mesh = mesh.translate([0.5 0.5 0.5]);
    bnd_mesh = mesh.bnd;
    
    % FEM Spaces
    Vcurl = fem(mesh,'NED');
    V = fem(mesh,'P1');
    
    % FEM Spaces with Dirichlet boundary conditions
    Vcurl0 = dirichlet(Vcurl,mesh.bnd);
    V0 = dirichlet(V,mesh.bnd);
    
    % Integration domain
    qud_order = 4;
    omega = dom(mesh,qud_order);
    
    % Assembling Galerkin matrix for the gauged formulation
    nu = 1; % Reluctivity
    curlcurlmat = nu * integral(omega,Vcurl0.curl,Vcurl0.curl);
    mixmat = integral(omega,Vcurl0,V0.grad);
    N_Vcurl0 = size(curlcurlmat,1);
    N_V0 = size(mixmat,2);
    zeromat = zeros(N_V0,N_V0);
    sysmat = [curlcurlmat mixmat; mixmat' zeromat];
    
    % Vector potential
    A = @(x,y,z) [cos(pi*x).*sin(pi*y).*sin(pi*z), -0.5*sin(pi*x).*cos(pi*y).*sin(pi*z), -0.5*sin(pi*x).*sin(pi*y).*cos(pi*z)];
    % Source 
    J = @(x,y,z) 3*pi^2*A(x,y,z);
    
    % Want to write an integral like integral(omega,Vcurl0,J)
    % Performing manual integration for pre_rhs_top
    
    % Getting the quadrature nodes and weights
    [X,W,elt2qud] = omega.qud;
    uqmat = Vcurl0.uqm(omega);
    JX = J(X(:,1),X(:,2),X(:,3));
    
    % RHS \int_{\Omega} J\cdot A' dx
    pre_rhs_top = sum(W.*(JX(:,1).*uqmat{1} + JX(:,2).*uqmat{2} + JX(:,3).*uqmat{3}),1)';
    N_temp = size(pre_rhs_top,1);
    rhs_top = eye(N_temp,N_temp)*pre_rhs_top;
    rhs_bot = zeros(N_V0,1);
    % RHS vector padded with zeros for the mixed formulation
    rhs_vec = [rhs_top; rhs_bot];
    
    % Solving the linear system
    sol = sysmat\rhs_vec;
    
    % Extracting the vector potential solution and the Lagrange multiplier
    sol_vecpot = sol(1:N_Vcurl0);
    sol_psi = sol(N_Vcurl0+1:end);
    
    % Projection to full fem spaces without Dirichlet BC
    P_Curl_Curl0 = elimination(Vcurl,mesh.bnd);
    sol_vecpot_full = P_Curl_Curl0 * sol_vecpot;
    P_V_V0 = elimination(V,mesh.bnd);
    sol_psi_full = P_V_V0 * sol_psi;
    
    full_dofs = Vcurl.dof;
    
    uqmat_full = Vcurl.uqm(omega);
    
    % Solution vector potential evaluated at quadrature nodes
    vecpot_sol_X = [sum(sol_vecpot_full'.*uqmat_full{1},2), sum(sol_vecpot_full'.*uqmat_full{2},2), sum(sol_vecpot_full'.*uqmat_full{3},2)];
    vecpot_X = A(X(:,1),X(:,2),X(:,3));
    
    vec_err_X = vecpot_sol_X - vecpot_X;
    err_X = dot(vec_err_X,vec_err_X,2);
    l2err = sum(W.*err_X,1);

    l2norm_sol = sum( dot(vecpot_sol_X,vecpot_sol_X,2) ,1);

    l2norm = sum( dot(vecpot_X,vecpot_X,2) ,1);

end