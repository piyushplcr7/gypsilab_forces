% Neumann problem solver DFK formulation

% TNA is passed as a vector of values of the neumann trace at qud pts of 
% Gamma

% NED is the space in which the dirichlet solution is sought

function [TDA] = Neumann_DFK(TNA,NED,method,order)
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

    % Ortho 2

% kernel = cell(3,1);
% kernel{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
% kernel{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
% kernel{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);
% 
% ortho2 = integral(Gamma,Gamma,NED,kernel,P1);
% ortho2 = ortho2 + regularize(Gamma,Gamma,NED,'grady[1/r]',P1);
% ortho2 = -ortho2/(4*pi);
% ortho = ortho2;
    
    NP1 = P1.ndof;
    Nned = NED.ndof;

    blockmat = [-Nmat ortho zeros(Nned,1);...
                ortho' zeros(NP1,NP1) B;...
                zeros(1,Nned) B' 0];

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
        rhs = [(-Cmat' + 0.5*M)*TNA_coeffs; zeros(NP1,1); 0];

        sol = blockmat\rhs;

        TDA = sol(1:Nned);

        disp('Error in discrete equation');
        norm(Nmat * TDA+ (-Cmat' + 0.5*M)*TNA_coeffs )
     


end