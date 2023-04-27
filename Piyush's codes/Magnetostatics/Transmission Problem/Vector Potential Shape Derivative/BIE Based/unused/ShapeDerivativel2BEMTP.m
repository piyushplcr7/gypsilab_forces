function val = ShapeDerivativel2BEMTP(bndmesh,g,Vel,omega_src,J)
    % 3rd term
    RWG = fem(bndmesh,'RWG');
    Gamma = dom(bndmesh,3);
    % Evaluating nxg at quadrature points where g lies in NED so nxg lies
    % in RWG
    nxgvals = reconstruct(g,Gamma,RWG);

    % Getting quadrature points for Gamma and omega_src
    [X,WX] = Gamma.qud;
    [Y,WY] = omega_src.qud;

    % Get tensor product quadrature rule
    NX = size(X,1);
    NY = size(Y,1);

    XX = repmat(X,NY,1); WWX = repmat(WX,NY,1); 
    nxgvals_XX = repmat(nxgvals,NY,1);
    YY = repelem(Y,NX,1); WWY = repelem(WY,NX,1);

    W = WWX .* WWY;

    % Kernel
    VelXX = Vel(XX);
    kernel = 3/(4*pi) * (XX-YY).*dot(XX-YY,VelXX,2)./vecnorm(XX-YY,2,2).^5 - 1/(4*pi)*VelXX./vecnorm(XX-YY,2,2).^3;
    JYY = J(YY(:,1),YY(:,2),YY(:,3));

    val = sum(W.*dot(cross(kernel,JYY,2),nxgvals_XX,2),1);
end