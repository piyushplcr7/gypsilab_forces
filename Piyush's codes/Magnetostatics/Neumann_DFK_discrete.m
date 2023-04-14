% Neumann problem solver DFK formulation

% TNA is passed as a vector of values of the neumann trace at qud pts of 
% Gamma

% NED is the space in which the dirichlet solution is sought

function [TDA] = Neumann_DFK_discrete(TNA,NED,method,order)
    bndmesh = NED.msh;
    Gamma = dom(bndmesh,order);

    % BEM spaces
    % curl conforming -> Dirichlet trace (NED)

    P1 = fem(bndmesh,'P1');

    % Kernel of the surface curl operator
    gradP1 = grad(P1); 

    % Div conforming space 
    RWG = fem(bndmesh,'RWG');

    %% LHS
    % Galerkin matrices

    % For operator N
    Nmat = -single_layer(Gamma,RWG.div,RWG.div);

    % Vector to enforce zero mean for P1 functions
    B = integral(Gamma,P1);

    % Matrix to enforce the orthogonal complement constraint
    ortho = single_layer(Gamma,NED,gradP1); % Uses my implementation.

    %% RHS

    switch method
        case "div"
        % For operator B, we can simply use C'
        % Operator C
        Cmat = double_layer_magnetostatics(Gamma,RWG,RWG);
        M = mass_matrix(Gamma,NED,RWG);
        TNA_coeffs = proj(TNA,Gamma,RWG);

        case "div0"
        DIV0 = nxgrad(P1);
        % For operator B, we can simply use C'
        % Operator C
        Cmat = double_layer_magnetostatics(Gamma,DIV0,RWG);
        M = mass_matrix(Gamma,NED,DIV0);
        TNA_coeffs = proj(TNA,Gamma,DIV0);
    end
        Nned = NED.ndof;
        %Nmat = (Nmat+Nmat')/2;
        % Solving in the orthogonal complement at discrete level
        [U,S,V] = svd(Nmat);
        nnzero = sum(diag(S) > 1e-14);
        % Determining the kernel
        kernel = V(:,nnzero+1:end);
        (Nmat*kernel);

        % Inner product matrix
        innerpdtmat = single_layer(Gamma,RWG.nx,RWG.nx) + single_layer(Gamma,RWG.div,RWG.div);
        %innerpdtmat = eye(Nned,Nned);

        % Augmented system
        lhsmat = [-Nmat innerpdtmat * kernel;
                   kernel' * innerpdtmat zeros(size(kernel,2),size(kernel,2))];

        rhsvec = [(-Cmat' + 0.5*M)*TNA_coeffs;
                   zeros(size(kernel,2),1)];

        soln = lhsmat\rhsvec;
        
        TDA = soln(1:Nned);

        % Printing out the part of solution that is expected to be zero
        %norm(soln(Nned+1:end))
        % Checking if the equation of concern is satisfied
        disp('Error in discrete equation');
        norm(Nmat * TDA+ (-Cmat' + 0.5*M)*TNA_coeffs )
        disp('Size of auxiliary solution (expected to be 0)');
        norm(soln(Nned+1:end))
        

end