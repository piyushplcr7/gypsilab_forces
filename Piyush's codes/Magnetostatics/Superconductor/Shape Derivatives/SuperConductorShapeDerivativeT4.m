function val = SuperConductorShapeDerivativeT4(bndmesh,TnA,DVel,omega_src,J)
    % 3rd term
    RWG = fem(bndmesh,'RWG');
    Gamma = dom(bndmesh,3);
    % Evaluating the Neumann Trace at the quadrature points
    Psi = reconstruct(TnA,Gamma,RWG);

    % Getting quadrature points for Gamma and omega_src
    [X,WX] = Gamma.qud;
   %  [Y,WY] = omega_src.qud;
   % 
   % % Get tensor product quadrature rule
   %  NX = size(X,1);
   %  NY = size(Y,1);
   % 
   %  XX = repmat(X,NY,1); WWX = repmat(WX,NY,1); Psi_XX = repmat(Psi,NY,1);
   %  YY = repelem(Y,NX,1); WWY = repelem(WY,NX,1);
   % 
   %  W = WWX .* WWY;
   % 
   %  % Evaluating rows of DVel for all quadrature points
   %  DVel1 = DVel{1}(XX);
   %  DVel2 = DVel{2}(XX);
   %  DVel3 = DVel{3}(XX);
   % 
   %  % Evaluating DVel * Psi
   %  DVel_Psi_XX = [dot(DVel1,Psi_XX,2) dot(DVel2,Psi_XX,2) dot(DVel3,Psi_XX,2)];
   % 
   %  % Kernel gradx G
   %  Gxy = @(X,Y) 1/(4*pi) * 1./vecnorm(X-Y,2,2);
   %  val = sum(W.* Gxy(XX,YY).*dot(DVel_Psi_XX,J(YY(:,1),YY(:,2),YY(:,3)),2),1 );

    % Computing using Gauss Quadrature
    data = omega_src.msh.col;
    R = data(1);
    r = data(2);
    AJ = computeVecpotTorus(1,R,r,X);
    val = sum(WX.* dot( [dot(DVel{1}(X),Psi,2) dot(DVel{2}(X),Psi,2) dot(DVel{3}(X),Psi,2)] , AJ , 2),1);
    
end