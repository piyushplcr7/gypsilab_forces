% This function implements the shape derivative assuming a constant
% magnetization

% Inputs
%
% TdA       : Coefficients for the exterior Dirichlet trace, in NED space
% TnA       : Coefficients for the exterior Neumann trace, in RWG

function sd = PermanentMagnetShapeDerivativeBIEVP(Gamma,TnA,TdA,J,omega_src,mu,mu0,M,Vel,DVel)
    %% BEM Spaces
    % Source mesh
    % Magnet boundary mesh
    bndmesh = Gamma.msh;
    [X,~] = Gamma.qud;
    normals = Gamma.qudNrm;
     % BEM spaces
    % curl conforming -> Dirichlet trace
    NED = fem(bndmesh,'NED'); 
    % Div conforming space 
    RWG = fem(bndmesh,'RWG');

    %% Remaining terms + partial derivative of l_M
    Nelt = bndmesh.nelt;
    [ii,jj] = meshgrid(1:Nelt,1:Nelt);

    Mvals = M(X);
    Mxn = cross(Mvals,normals,2);
    Mxn_coeffs = proj(Mxn,Gamma,RWG);
    
    kernelA1 = @(x,y,z) dot(z,Vel(x) - Vel(y), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    kernelA2 = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
    kernelC1 = @(x,y,z) 1/(4*pi) * z./vecnorm(z,2,2).^3;
    kernelC3 = @(x,y,z) -3/(4*pi) * z .* dot(z,Vel(y)-Vel(x),2)./vecnorm(z,2,2).^5 + 1/(4*pi)*(Vel(y)-Vel(x))./vecnorm(z,2,2).^3;
    kernelN = kernelA1;

    DVelRWG = RWG;
    DVelRWG.opr = 'Dvel[psi]';

    euler = parcluster('local');
    euler.NumWorkers = 5;
    saveProfile(euler);

    pool = euler.parpool(5);

    spmd
        if spmdIndex==1
            % Red term in implementation.xopp, pg 67
            % z := y-x, grad_x G = z/norm(z)^3
            A1mat = panel_assembly(bndmesh,kernelA1,RWG,RWG,ii(:),jj(:));
            red_remaining = mu/2 * Mxn_coeffs' * A1mat * Mxn_coeffs;
            red_l_M = mu * TnA' * A1mat * Mxn_coeffs;
            A1 = TnA' * A1mat * TnA;

        elseif spmdIndex==2
            % Blue terms in implementation.xopp, pg 67
            A2mat = panel_assembly_shape_derivative(bndmesh,kernelA2,DVelRWG,RWG,ii(:),jj(:),Vel,DVel);
            blue_remaining = mu * Mxn_coeffs' * A2mat * Mxn_coeffs;
            blue_l_M = mu * (Mxn_coeffs' * A2mat * TnA + TnA' * A2mat * Mxn_coeffs);
            A2 = 2* TnA' * A2mat * TnA;

        elseif spmdIndex==3
            % Derivative of double layer in l_M, borrowed from TP
            % Partial derivative of b_C
            % z := y-x
            C1mat = panel_assembly_shape_derivative(bndmesh,kernelC1,DVelRWG,RWG,ii(:),jj(:),Vel,DVel);
            MC1 = -mu0* Mxn_coeffs' * C1mat * TdA;
            MC2 = -mu0* TdA' * C1mat * Mxn_coeffs; 
            C1 = TnA' * C1mat * TdA;
            C2 = TdA' * C1mat * TnA; 

        elseif spmdIndex==4
            % C3 (Is this way of evaluation okay?), z:= y-x
            C3mat = panel_assembly(bndmesh,kernelC3,RWG,RWG,ii(:),jj(:));
            MC3 = -mu0* Mxn_coeffs' * C3mat * TdA;
            C3 = TnA' * C3mat * TdA;

        elseif spmdIndex==5
            % Partial derivative of b_N
            Nmat = panel_assembly_shape_derivative(bndmesh,kernelN,RWG.div,RWG.div,ii(:),jj(:),Vel,DVel);
            N = -TdA' * Nmat * TdA;

        end

    end

    %% Shape derivative of linear forms from transmission problem

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
    nxgvals = cross(normals,gvals,2);
    nxgvals_XX = repmat(nxgvals,NY,1);
    kernell22 = 3/(4*pi) * (XX-YY).*dot(XX-YY,VelXX,2)./vecnorm(XX-YY,2,2).^5 - 1/(4*pi)*VelXX./vecnorm(XX-YY,2,2).^3;
    l22 = sum(W.*dot(cross(kernell22,JYY,2),nxgvals_XX,2),1);

    % Kernel gradx G
    gradxG = 1/(4*pi) * (YY-XX)./vecnorm(XX-YY,2,2).^3;
    % DVel {nxg}
    DVelnxgvalsXX = [dot(DVel1XX,nxgvals_XX,2) dot(DVel2XX,nxgvals_XX,2) dot(DVel3XX,nxgvals_XX,2) ];
    l21 = sum(W.* dot(cross(gradxG,JYY,2),DVelnxgvalsXX,2),1);

    % SD from transmission problem
    sd = 1/(2*mu0) * ( (1+mu/mu0) * (A1{1}+A2{2})...
                        -4 * (C1{3}+C2{3}+C3{4})...
                        +(1+mu0/mu) * (N{5}) )...
                        -l1 + (l21+l22); 
    % SD additional terms for magnet
    sd = sd + red_remaining{1} + red_l_M{1} + blue_remaining{2} + blue_l_M{2} + MC1{3} + MC2{3} + MC3{4};
    pool.delete();
end