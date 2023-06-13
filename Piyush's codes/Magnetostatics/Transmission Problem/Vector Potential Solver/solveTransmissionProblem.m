% Solve Transmission Problem. Returns the exterior traces 

function [Psi,g] = solveTransmissionProblem(bndmesh,J,omega_src,mu_i,mu_e)
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
    % LHS Matrix
%     blockopr = [(1+mu_i/mu_e)*Amat -2*Cmat vec;
%                 2*Bmat -(1+mu_e/mu_i)*Nmat zeros(Nned,1);
%                 vec' zeros(1,Nned) 0];

    blockopr = [(1+mu_i/mu_e)*Amat -2*Cmat zeros(NP1,NP1) vec zeros(NP1,1);
                2*Bmat -(1+mu_e/mu_i)*Nmat ortho zeros(Nned,1) zeros(Nned,1);
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
    TNAJ_DIV_coeffs = proj(TNAJ,Gamma,DIV);
%     TNAJ_DIV0_coeffs = proj(TNAJ,Gamma,DIV0);

    M_div0_ned = mass_matrix(Gamma,DIV0,NED);
    M_div_ned = mass_matrix(Gamma,DIV,NED);

    % mu_e <Td^+ N(J), zeta>
    rhs1 = mu_e * M_div0_ned * TDAJ_NED_coeffs;
    % mu_e <Tn^+ N(J), u >
%     rhs2 = mu_e * M_div0_ned' * TNAJ_DIV0_coeffs;
    rhs2 = mu_e * M_div_ned' * TNAJ_DIV_coeffs;

    rhs = [rhs1; rhs2; zeros(NP1,1); 0 ;0];
    sol = blockopr\rhs;

    Psi = sol(1:NP1);
    g = sol(NP1+1:NP1+Nned);

end