% TnA and TdA are the coeffs of the exterior traces of the vector potential
% Psi and g respectively. TnA is expected to be in RWG space.

function sd = SdBemTPVP(bndmesh,TnA,TdA,J,omega_src,Vel,DVel)
    % BEM Spaces
    P1 = fem(bndmesh,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0 = nxgrad(P1); 
    RWG = fem(bndmesh,'RWG');
    NED = fem(bndmesh,'NED');
    Gamma = dom(bndmesh,3);

    Nelt = bndmesh.nelt;
    [ii,jj] = meshgrid(1:Nelt,1:Nelt);

    % partial derivative of b_A 
    kernelA1 = @(x,y,z) dot(z,Vel(x) - Vel(y), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    A1mat = panel_assembly(bndmesh,kernelA1,RWG,RWG,ii(:),jj(:));
    A1 = TnA' * A1mat * TnA;
    kernelA2 = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
    %
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
    kernelN = @(x,y,z)


end