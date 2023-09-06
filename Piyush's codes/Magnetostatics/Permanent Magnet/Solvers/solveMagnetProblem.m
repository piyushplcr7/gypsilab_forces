% BIE solver for a permanent magnet

function [TnA,TdA] = solveMagnetProblem(mesh,bndmesh,J,omega_src,mu,mu0,M)
    %% BEM Spaces
    Gamma = dom(bndmesh,3);
    [X,~] = Gamma.qud;
    % BEM spaces
    % curl conforming -> Dirichlet trace
    NED = fem(bndmesh,'NED'); 
    P1 = fem(bndmesh,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0 = nxgrad(P1); 
    % Kernel of the surface curl operator
    Ker_curl = grad(P1); 
    % Div conforming space 
    DIV = fem(bndmesh,'RWG');

    %% Galerkin matrices 
    % For operator A
    Amat = single_layer(Gamma,DIV0,DIV0);
    % For operator C
    Cmat = double_layer_magnetostatics(Gamma,DIV0,DIV);
    % For operator B, we can simply use C'
    Bmat = Cmat';
    % For operator N
    Nmat = -single_layer(Gamma,DIV.div,DIV.div);

    NP1 = size(Amat,1); % Equal to NP1
    Nned = size(Nmat,1);

    % vector to enforce zero mean for P1 functions
    vec = integral(Gamma,P1);

    ortho = single_layer(Gamma,NED,Ker_curl); % Uses my implementation.

    %% Modified Linear System

    blockopr = [(1+mu/mu0)*Amat -2*Cmat zeros(NP1,NP1) vec zeros(NP1,1);
                2*Bmat -(1+mu0/mu)*Nmat ortho zeros(Nned,1) zeros(Nned,1);
                zeros(NP1,NP1) ortho' zeros(NP1,NP1) zeros(NP1,1) vec;
                vec' zeros(1,Nned) zeros(1,NP1) 0 0;
                zeros(1,NP1) zeros(1,Nned) vec' 0 0];

    % Computing the RHS
    % Computing the fields on the points X
    AJ = compute_vecpot(J,omega_src,X);
    curlAJ = compute_vecpot_curl(J,omega_src,X);

    % Taking traces of the fields at the points X
    normals = Gamma.qudNrm;
    TDAJ = AJ - dot(AJ,normals,2).*normals;
    TNAJ = cross(curlAJ,normals,2);

    % Projecting the Dirichlet trace to NED space
    TDAJ_NED_coeffs = proj(TDAJ,Gamma,NED);

    % Projecting the Neumann trace to DIV space
    %TNA_DIV_coeffs = proj(TNA,Gamma,DIV);
    TNAJ_DIV0_coeffs = proj(TNAJ,Gamma,DIV0);

    M_div0_ned = mass_matrix(Gamma,DIV0,NED);

    % Adding the magnetization induced RHS terms
    % Projecting Mxn to RWG
    Mvals = M(X);
    Mxn = cross(Mvals,normals,2);
    Mxn_coeffs = proj(Mxn,Gamma,DIV0);

    % Computing the interacting bnd-vol integrals
    KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
    KDL = @(x,y,z)1/(4*pi)*z./vecnorm(z,2,2).^3;
    NEDVol = fem(mesh,'NED');
    Nelt_bnd = bndmesh.nelt;
    Nelt_vol = mesh.nelt;
    [II,JJ] = meshgrid(1:Nelt_bnd,1:Nelt_vol);
    intmat1 = panel_assembly_cross_modif(bndmesh,mesh,KV,NEDVol.curl,RWG,II(:),JJ(:));
    intmat2 = panel_assembly_cross_modif(bndmesh,mesh,KDL,NEDVol.curl,RWG,II(:),JJ(:));
    % Projecting curlM to NEDVol.curl by projecting M to NEDVol
    omega = dom(mesh,4);
    [Xvol,~] = omega.qud;
    Mvol = M(Xvol);
    Mvol_coeffs = proj(Mvol,mesh,NEDVol);

    % mu_e <Td^+ N(J), zeta>
    rhs1 = mu0 * M_div0_ned * TDAJ_NED_coeffs - mu * Amat * Mxn_coeffs...
            - mu * intmat1 * Mvol_coeffs;
    % mu_e <Tn^+ N(J), u >
    rhs2 = mu0 * M_div0_ned' * TNAJ_DIV0_coeffs...
            + mu0/2 * M_div0_ned' * Mxn_coeffs...
            - mu0 * Bmat * Mxn_coeffs - mu0 * intmat2 * Mvol_coeffs;

    rhs = [rhs1; rhs2; zeros(NP1,1); 0 ;0];
    sol = blockopr\rhs;

    TnA = sol(1:NP1);
    TdA = sol(NP1+1:NP1+Nned);
end