% Magnetostatics BEM implementation test

addpath(genpath("../../"));

clear; clc;

for i = 3:10
    N = 2^i;
    disp(N);
    %N = 50;
    
    L = 5*[1 1 1];
    
    qud_order_vol = 3;
    qud_order_bnd = 3;
    
    % Solution domain
    mesh = mshCube(N,L);
    bndmesh = mesh.bnd;
    bndmesh = bndmesh.translate([10 0 0]);
    
    % Source parameters
    % Geometry
    N_src = N;
    R0 = 2;
    r0 = .5;
    mesh_src = mshTorus(N_src,R0,r0);
    
    omega_src = dom(mesh_src,qud_order_vol);
    Gamma = dom(bndmesh,qud_order_bnd);
    
    % Current
    interior_torus = @(x,y,z) sqrt( (x-x*R0./(sqrt(x.^2+y.^2))).^2 + (y-y*R0./(sqrt(x.^2+y.^2))).^2 + z.^2 )<r0;
    torus_tangent = @(x,y,z) interior_torus(x,y,z).*[-y x z*0]./(sqrt(x.^2+y.^2));
    J0 = 1;
    J = @(x,y,z) J0*torus_tangent(x,y,z); 
    
    % Combined mesh visualized
    mesh_combined = union(mesh_src,bndmesh);
    plot(mesh_combined);
    
    %%
    
    Vtrial = fem(bndmesh,'NED');

    % RWG is not div0 space
    %Vtest = fem(bndmesh,'RWG');

    % Getting the div0 space as the surface curl or nxgrad of P1
    P1 = fem(bndmesh,'P1');
    Vtest = nxgrad(P1);
    
    % Need to evaluate the Dirichlet and Neumann trace of the vector potential
    % generated by a current inside the torus (mesh_src)
    Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); 
    GradGy = cell(3,1);
    GradGy{1} = @(X,Y)1/4/pi*femGreenKernel(X,Y,'grady[1/r]1',0);
    GradGy{2} = @(X,Y)1/4/pi*femGreenKernel(X,Y,'grady[1/r]2',0);
    GradGy{3} = @(X,Y)1/4/pi*femGreenKernel(X,Y,'grady[1/r]3',0);

    %%
    % BEM Galerkin matrices
    %A = single_layer(bndmesh,qud_order_bnd);
    C = double_layer(bndmesh,qud_order_bnd);
    M = integral(Gamma,Vtest,Vtrial);

    % Computing the Single Layer matrix A
    A = integral(Gamma,Gamma,Vtest,Gxy,Vtest)/(4*pi);
    A = A + 1/(4*pi)*regularize(Gamma,Gamma,Vtest,'[1/r]',Vtest);

    % This matrix is singular, regularize it
    P0 = fem(bndmesh,'P0');
    V = 1/4/pi*integral(Gamma,Gamma,P0,Gxy,P0);
    V = V + 1/4/pi*regularize(Gamma,Gamma,P0,'[1/r]',P0);

    Nv = size(V,1);
    M_beta = integral(Gamma,P0,P0);

    eq_density = V\(M_beta * ones(Nv,1));

    M_mixed = integral(Gamma,P1,P0);

    alpha = 0.01;

    lhsmat = A + alpha * M_mixed * eq_density * eq_density' * M_mixed';


    %%
    [X_src,W_src,elt2qud_src] = omega_src.qud;
    [X,W,elt2qud] = Gamma.qud;
    
    J_X_src = J(X_src(:,1),X_src(:,2),X_src(:,3));
    % Need to project the Dirichlet trace to Vtrial as well!!!!
    
    N_qud_pts = size(X,1);
    N_qud_pts_src = size(X_src,1);
    
    integrand_without_green_fn = W_src.*J_X_src;
    
    % Repeating the source quadrature points N_qud_pts times
    X_src_extended = repmat(X_src,N_qud_pts,1);
    X_extended = repelem(X,N_qud_pts_src,1);
    
    green_fn = Gxy(X_src_extended,X_extended);
    green_fn = reshape(green_fn,[N_qud_pts_src, N_qud_pts]);
    
    % Vector potential evaluated at X
    integral1 = sum(integrand_without_green_fn(:,1).*green_fn,1)';
    integral2 = sum(integrand_without_green_fn(:,2).*green_fn,1)';
    integral3 = sum(integrand_without_green_fn(:,3).*green_fn,1)';
    
    % Taking the trace of the vector potential on bndmesh
    normals = bndmesh.nrm;
    normals_X = repelem(normals,3,1);
    vecpot_X = [integral1 integral2 integral3];
    TD_vecpot_X = vecpot_X - dot(normals_X,vecpot_X,2).*normals_X;
    
    integral1 = TD_vecpot_X(:,1);
    integral2 = TD_vecpot_X(:,2);
    integral3 = TD_vecpot_X(:,3);
    
    % Visualizing the Dirichlet trace of the vector potential
    figure;
    quiver3(X(:,1),X(:,2),X(:,3),integral1,integral2,integral3);
    
    % Project the computed vector potential to Vtrial and get coeffs TD
    M_tr_tr = integral(Gamma,Vtrial,Vtrial);
    % RHS for projection
    uqmat = Vtrial.uqm(Gamma);
    rhs = sum(W.*(uqmat{1}.*integral1+uqmat{2}.*integral2+uqmat{3}.*integral3),1)';
    TD = M_tr_tr\rhs;
    
    % Getting the Neumann Trace
    % Evaluating the components of gradient of greens function
    grad_green_fn_1 = GradGy{1}(X_src_extended,X_extended);
    grad_green_fn_2 = GradGy{2}(X_src_extended,X_extended);
    grad_green_fn_3 = GradGy{3}(X_src_extended,X_extended);
    
    % Reshaping
    grad_green_fn_1 = reshape(grad_green_fn_1,[N_qud_pts_src,N_qud_pts]);
    grad_green_fn_2 = reshape(grad_green_fn_2,[N_qud_pts_src,N_qud_pts]);
    grad_green_fn_3 = reshape(grad_green_fn_3,[N_qud_pts_src,N_qud_pts]);
    
    % Computing the cross product in the integrand
    blabla1 = J_X_src(:,2).*grad_green_fn_3 - J_X_src(:,3).*grad_green_fn_2;
    blabla2 = J_X_src(:,3).*grad_green_fn_1 - J_X_src(:,1).*grad_green_fn_3;
    blabla3 = J_X_src(:,1).*grad_green_fn_2 - J_X_src(:,2).*grad_green_fn_1;
    
    % Computing the components of the integral, in the end it contains the curl
    % of the vector potential at the bndmesh qud points
    blabla1 = sum(W_src.*blabla1,1)';
    blabla2 = sum(W_src.*blabla2,1)';
    blabla3 = sum(W_src.*blabla3,1)';
    
    curl_vecpot_X = [blabla1 blabla2 blabla3];
    
    % Computing n cross curl_vecpot_X
    TN_vecpot_X = cross(normals_X,curl_vecpot_X,2);
    
    blabla1 = TN_vecpot_X(:,1);
    blabla2 = TN_vecpot_X(:,2);
    blabla3 = TN_vecpot_X(:,3);
    
    % Projecting to space Vtest and getting the coefficients
    M_ts_ts = integral(Gamma,Vtest,Vtest);
    % RHS for projection
    uqmat_test = Vtest.uqm(Gamma);
    rhs_test = sum(W.*(uqmat_test{1}.*blabla1+uqmat_test{2}.*blabla2+uqmat_test{3}.*blabla3),1)';
    TN_known = M_ts_ts\rhs_test;
    
    % Visualizing the Neumann Trace of the vector potential
    figure;
    quiver3(X(:,1),X(:,2),X(:,3),blabla1,blabla2,blabla3);
    
    %% Solving for the Neumann Trace
    TN = lhsmat\((0.5*M+C)*TD);
    
    %% Comuting the errors
    
    err_vec = TN-TN_known;
    
    % L2 error
    l2err = err_vec'*M_ts_ts*err_vec
    
    % Energy norm error? (by A)
    energerr = err_vec'*A*err_vec
    close all;
end