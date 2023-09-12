% BIE solver for a permanent magnet. Returns the solution split into two
% parts

function [TnA,TdA] = solveTransProblemConstField(Gamma,mu,mu0,B0)
    %% BEM Spaces
    jumpNu = 1/mu0-1/mu;
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

    normals = Gamma.qudNrm;

    M_div0_ned = mass_matrix(Gamma,DIV0,NED);

    % Adding the B0 induced RHS terms
    % Projecting B0xn to RWG
    B0vals = B0(X);
    B0xn = cross(B0vals,normals,2);
    B0xn_coeffs = projSpecial(B0xn,Gamma,DIV0);

    % RHS1 from B0xn
    rhsM1 = - mu * jumpNu * Amat * B0xn_coeffs;
    % RHS2 from B0xn
    rhsM2 = mu0/2 * jumpNu * M_div0_ned' * B0xn_coeffs...
            - mu0 * jumpNu * Bmat * B0xn_coeffs;

    rhsM = [rhsM1; rhsM2; zeros(NP1,1); 0 ;0];

    solM = blockopr\rhsM;

    TnA = solM(1:NP1);
    TdA = solM(NP1+1:NP1+Nned);
end