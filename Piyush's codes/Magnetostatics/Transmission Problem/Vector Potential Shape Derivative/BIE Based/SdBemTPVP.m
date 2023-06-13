% TnA and TdA are the coeffs of the exterior traces of the vector potential
% Psi and g respectively. TnA is expected to be in RWG space.

function sd = SdBemTPVP(bndmesh,TnA,TdA,J,omega_src,Vel,DVel,mu0,mu)
    order = 3;
    % BEM Spaces
    P1 = fem(bndmesh,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0 = nxgrad(P1); 
    RWG = fem(bndmesh,'RWG');
    NED = fem(bndmesh,'NED');
    Gamma = dom(bndmesh,order);

    Nelt = bndmesh.nelt;
    [ii,jj] = meshgrid(1:Nelt,1:Nelt);

    % partial derivative of b_A 
    kernelA1 = @(x,y,z) dot(z,Vel(x) - Vel(y), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    A1mat = panel_assembly(bndmesh,kernelA1,RWG,RWG,ii(:),jj(:));
    A1 = TnA' * A1mat * TnA;
    kernelA2 = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
    %g
    DVelRWG = RWG;
    DVelRWG.opr = 'Dvel[psi]';
    A2mat = panel_assembly_shape_derivative(bndmesh,kernelA2,DVelRWG,RWG,ii(:),jj(:),Vel,DVel);
    A2 = 2* TnA' * A2mat * TnA;

    % Partial derivative of b_C
    % z := y-x
    kernelC1 = @(x,y,z) 1/(4*pi) * z./vecnorm(z,2,2).^3;
    C1mat = panel_assembly_shape_derivative(bndmesh,kernelC1,DVelRWG,RWG,ii(:),jj(:),Vel,DVel);
    C1 = TnA' * C1mat * TdA;
    C2 = TdA' * C1mat * TnA; % = C1?
    % C3 (Is this way of evaluation okay?), z:= y-x
    kernelC3 = @(x,y,z) -3/(4*pi) * z .* dot(z,Vel(y)-Vel(x),2)./vecnorm(z,2,2).^5 + 1/(4*pi)*(Vel(y)-Vel(x))./vecnorm(z,2,2).^3;
    C3mat = panel_assembly(bndmesh,kernelC3,RWG,RWG,ii(:),jj(:));
    C3 = TnA' * C3mat * TdA;

    % Partial derivative of b_N
    kernelN = kernelA1;
    Nmat = panel_assembly_shape_derivative(bndmesh,kernelN,RWG.div,RWG.div,ii(:),jj(:),Vel,DVel);
    N = -TdA' * Nmat * TdA;

    % Partial derivative of l1, computed using superconductor shape
    % derivative implementation
    l1 = SuperConductorShapeDerivativeT3(bndmesh,TnA,Vel,omega_src,J)...
        + SuperConductorShapeDerivativeT4(bndmesh,TnA,DVel,omega_src,J);

    % Partial derivaive of l2
    % Getting quadrature points for Gamma and omega_src
    [X,WX] = Gamma.qud;
    [Y,WY] = omega_src.qud;
    % Get tensor product quadrature rule
    NX = size(X,1);
    NY = size(Y,1);
    XX = repmat(X,NY,1); WWX = repmat(WX,NY,1); 
    YY = repelem(Y,NX,1); WWY = repelem(WY,NX,1);
    W = WWX .* WWY;
    VelXX = Vel(XX);
    DVel1XX = DVel{1}(XX);
    DVel2XX = DVel{2}(XX);
    DVel3XX = DVel{3}(XX);
    JYY = J(YY(:,1),YY(:,2),YY(:,3));

    gvals = reconstruct(TdA,Gamma,NED);
    normals = Gamma.qudNrm;
    nxgvals = cross(normals,gvals,2);
    nxgvals_XX = repmat(nxgvals,NY,1);
    kernell22 = 3/(4*pi) * (XX-YY).*dot(XX-YY,VelXX,2)./vecnorm(XX-YY,2,2).^5 - 1/(4*pi)*VelXX./vecnorm(XX-YY,2,2).^3;
    l22 = sum(W.*dot(cross(kernell22,JYY,2),nxgvals_XX,2),1);

    % Kernel gradx G
    gradxG = 1/(4*pi) * (YY-XX)./vecnorm(XX-YY,2,2).^3;
    % DVel {nxg}
    DVelnxgvalsXX = [dot(DVel1XX,nxgvals_XX,2) dot(DVel2XX,nxgvals_XX,2) dot(DVel3XX,nxgvals_XX,2) ];
    l21 = sum(W.* dot(cross(gradxG,JYY,2),DVelnxgvalsXX,2),1);

    sd = 1/(2*mu0) * ( (1+mu/mu0) * (A1+A2)...
                        -4 * (C1+C2+C3)...
                        +(1+mu0/mu) * (N) )...
                        -l1 + (l21+l22);

end