% Debug solveTpPMCFSP

function [psi_i,g_i,psi_e] = solveTpPMCFSP_debug(bndmesh_i,bndmesh_e,mu0,M)
    % BEM Spaces
    P1_i = fem(bndmesh_i,'P1');
    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');

    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    normals_i = Gamma_i.qudNrm;

    % Computing M . n at Gamma_i
    [X_i,~] = Gamma_i.qud;
    Mvals = M(X_i);
    Mdotn = dot(Mvals,normals_i,2);
    Mdotncoeffs = proj(Mdotn,Gamma_i,P0_i);

    %% Operator matrices

    Vii = single_layer(Gamma_i,P0_i,P0_i);
    Vee = single_layer(Gamma_e,P0_e,P0_e);

    Kii = double_layer_laplace(Gamma_i,P0_i,P1_i);
    Wii = single_layer(Gamma_i,P1_i.nxgrad,P1_i.nxgrad);
    
    % Cross matrices
    Vei = single_layer_cross(Gamma_i,Gamma_e,P0_i,P0_e);
    Kie = double_layer_laplace_cross(Gamma_e,Gamma_i,P0_e,P1_i);

    %% Linear system

    blockopr = [2*Vii, 2*Kii, Vei;
                -2*Kii', 2*Wii, -Kie';
                Vei', Kie, Vee];

    rhs = [-Vii * Mdotncoeffs;
            -0.5 * mass_matrix(Gamma_i,P1_i,P0_i) * Mdotncoeffs + Kii' * Mdotncoeffs;
            zeros(P0_e.ndof,1)];

    sol = blockopr\rhs;

    psi_i = sol(1:P0_i.ndof);
    g_i = sol(P0_i.ndof+1:P0_i.ndof+P1_i.ndof);
    psi_e = sol(P0_i.ndof+P1_i.ndof+1:P0_i.ndof+P1_i.ndof+P0_e.ndof);

    %% Debug
    % Computing the trace of -psi_sl(M.n) at the outer boundary
    [X_e,~] = Gamma_e.qud;
    [Y_i,W_i] = Gamma_i.qud;

    NX = size(X_e,1);
    NY = size(Y_i,1);

    XX = repelem(X_e,NY,1);
    YY = repmat(Y_i,NX,1);
    WYY_i = repmat(W_i,NX,1);
    MdotnYY_i = repmat(Mdotn,NX,1);

    integrand = -1/4/pi./vecnorm(XX-YY,2,2).*MdotnYY_i;
    integrand = reshape(integrand,[NY,NX]);
    % Traces at the outer boundary of -psi_sl(M.n)
    g_w_vals = sum(W_i.*integrand,1)';
    P1_e = fem(bndmesh_e,'P1');
    g_w_coeffs = proj(g_w_vals,Gamma_e,P1_e);
    Kee = double_layer_laplace(Gamma_e,P0_e,P1_e);

    psi_v = Vee\(0.5*mass_matrix(Gamma_e,P0_e,P1_e)*g_w_coeffs + Kee * g_w_coeffs);


    % Traces of the representation formula for v on gamma_i
    [X_i,~] = Gamma_i.qud;
    [Y_e,W_e] = Gamma_e.qud;
    psi_v_vals = reconstruct(psi_v,Gamma_e,P0_e);
    normals_e = Gamma_e.qudNrm;

    NX = size(X_i,1);
    NY = size(Y_e,1);

    XX = repelem(X_i,NY,1);
    YY = repmat(Y_e,NX,1);
    normals_e_YY = repmat(normals_e,NX,1);

    Gxy  = reshape(1/4/pi./vecnorm(XX-YY,2,2),[NY,NX]);
    gradyGdotn = reshape(1/4/pi./vecnorm(XX-YY,2,2).^3.*dot(XX-YY,normals_e_YY,2),[NY,NX]);

    Tdpsislpsiv = sum(W_e.*Gxy.*psi_v_vals,1)';
    yo1coeffs = proj(Tdpsislpsiv,Gamma_i,P1_i);
    Tdpsidlgw = sum(W_e.*gradyGdotn.*g_w_vals,1)';
    yo2coeffs = proj(Tdpsidlgw,Gamma_i,P1_i);

    % Comparing the constructed solution with solved g_i
    % g_i = -V(Mdotn) + psi_sl(psi_v) - psi_dl(g_w)

    M01i = mass_matrix(Gamma_i,P0_i,P1_i);
    lhs = M01i * g_i;
    rhs = -Vii * Mdotncoeffs + M01i * yo1coeffs - M01i * yo2coeffs;
    plot(lhs)
    hold on
    plot(rhs)
    disp("what");


end