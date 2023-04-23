% Shape derivative computation of force (simplified)

function sd = ForceSdBemTP(bndmesh,TnA,g,J,omega_src,Vel)
    % BEM Spaces
    P1 = fem(bndmesh,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0 = nxgrad(P1); 
    RWG = fem(bndmesh,'RWG');
    NED = fem(bndmesh,'NED');
    Gamma = dom(bndmesh,3);

    % Evaluating the Neumann Trace at the quadrature points
    Psi = reconstruct(TnA,Gamma,DIV0);

    % Getting quadrature points for Gamma and omega_src
    [X,WX] = Gamma.qud;
    [Y,WY] = omega_src.qud;

    % Get tensor product quadrature rule
    NX = size(X,1);
    NY = size(Y,1);

    XX = repmat(X,NY,1); WWX = repmat(WX,NY,1); Psi_XX = repmat(Psi,NY,1);
    YY = repelem(Y,NX,1); WWY = repelem(WY,NX,1);

    W = WWX .* WWY;

    VelXX = Vel(XX);
    JYY = J(YY(:,1),YY(:,2),YY(:,3));
    
    % Kernel gradx G
    gradxG = @(X,Y) 1/(4*pi) * (Y-X)./vecnorm(X-Y,2,2).^3;
    dl1_ds = sum(W.* dot(gradxG(XX,YY),VelXX,2).*dot(Psi_XX,JYY,2),1 );


    %% 
    gvals = reconstruct(g,Gamma,NED);
    normals = Gamma.qudNrm;
    nxgvals = cross(normals,gvals,2);
    nxgvals_XX = repmat(nxgvals,NY,1);
    kernel = 3/(4*pi) * (XX-YY).*dot(XX-YY,VelXX,2)./vecnorm(XX-YY,2,2).^5 - 1/(4*pi)*VelXX./vecnorm(XX-YY,2,2).^3;
    dl2_ds = sum(W.*dot(cross(kernel,JYY,2),nxgvals_XX,2),1);

    sd = -dl1_ds + dl2_ds;
end