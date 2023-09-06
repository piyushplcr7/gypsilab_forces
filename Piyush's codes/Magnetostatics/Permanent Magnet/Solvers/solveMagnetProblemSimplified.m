% BIE solver for a permanent magnet. Returns the solution split into two
% parts

function [TnAJ,TdAJ,TnAM,TdAM] = solveMagnetProblemSimplified(Gamma,J,omega_src,mu,mu0,M)
    %% BEM Spaces
    bndmesh = Gamma.msh;
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
    %TNAJ_DIV0_coeffs = proj(TNAJ,Gamma,DIV0);
    TNAJ_DIV0_coeffs = projSpecial(TNAJ,Gamma,DIV0);

    M_div0_ned = mass_matrix(Gamma,DIV0,NED);

    % Adding the magnetization induced RHS terms
    % Projecting Mxn to RWG
    Mvals = M(X);
    Mxn = cross(Mvals,normals,2);
%     Mxn_coeffs = proj(Mxn,Gamma,DIV0);
    Mxn_coeffs = projSpecial(Mxn,Gamma,DIV0);

    % mu_e <Td^+ N(J), zeta>
    rhsJ1 = mu0 * M_div0_ned * TDAJ_NED_coeffs;
    % mu_e <Tn^+ N(J), u >
    rhsJ2 = mu0 * M_div0_ned' * TNAJ_DIV0_coeffs;

    % Magnetic RHS
    rhsM1 = - mu * Amat * Mxn_coeffs;
    % Magnetic RHS
    rhsM2 = mu0/2 * M_div0_ned' * Mxn_coeffs...
            - mu0 * Bmat * Mxn_coeffs;

    rhsJ = [rhsJ1; rhsJ2; zeros(NP1,1); 0 ;0];
    rhsM = [rhsM1; rhsM2; zeros(NP1,1); 0 ;0];

    solJ = blockopr\rhsJ;
    solM = blockopr\rhsM;

    TnAJ = solJ(1:NP1);
    TdAJ = solJ(NP1+1:NP1+Nned);

    TnAM = solM(1:NP1);
    TdAM = solM(NP1+1:NP1+Nned);
end