function val = SuperConductorShapeDerivativeT3(bndmesh,TnA,Vel,omega_src,J)
    % 3rd term
    RWG = fem(bndmesh,'RWG');
    Gamma = dom(bndmesh,3);
    % Evaluating the Neumann Trace at the quadrature points
    Psi = reconstruct(TnA,Gamma,RWG);

    % Getting quadrature points for Gamma and omega_src
    [X,WX] = Gamma.qud;
    % [Y,WY] = omega_src.qud;
    % 
    % % Get tensor product quadrature rule
    % NX = size(X,1);
    % NY = size(Y,1);
    % 
    % XX = repmat(X,NY,1); WWX = repmat(WX,NY,1); Psi_XX = repmat(Psi,NY,1);
    % YY = repelem(Y,NX,1); WWY = repelem(WY,NX,1);
    % 
    % W = WWX .* WWY;
    % 
    % % Kernel gradx G
    % gradxG = @(X,Y) 1/(4*pi) * (Y-X)./vecnorm(X-Y,2,2).^3;
    % val = sum(W.* dot(gradxG(XX,YY),Vel(XX),2).*dot(Psi_XX,J(YY(:,1),YY(:,2),YY(:,3)),2),1 );

    % Computing with Gaussian quadrature for the integration over source
    data = omega_src.msh.col;
    R = data(1);
    r = data(2);

    middle = computeGradxGJTTorus(1,R,r,X);
    val = sum(WX.*dot(Vel(X), [dot(middle{1}, Psi,2) dot(middle{2}, Psi,2) dot(middle{3}, Psi,2)] , 2 ), 1);
end