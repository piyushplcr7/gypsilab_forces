function sd = SdBemLMCFVP(bndmesh_i,bndmesh_e,Psi_i,g_i,Psi_e,Vel,DVel,mu0,mu,B0)
    jumpMuInv = 1/mu0-1/mu;
    % Integration domain
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    % BEM Spaces for Gamma_i
    P1_i = fem(bndmesh_i,'P1');
    DIV0_i = nxgrad(P1_i);
    RWG_i = fem(bndmesh_i,'RWG');

    % BEM Spaces for gamma_e
    P1_e = fem(bndmesh_e,'P1');
    DIV0_e = nxgrad(P1_e);
    RWG_e = fem(bndmesh_e,'RWG');

    Psi_i_vals = reconstruct(Psi_i,Gamma_i,RWG_i);
    Psi_e_vals = reconstruct(Psi_e,Gamma_e,RWG_e);
    nxgvals_i = reconstruct(g_i,Gamma_i,RWG_i);

    % Projecting B0xn to RWG_i space
    normals_i = Gamma_i.qudNrm;
    % B0xn at quadrature points
    B0xn = cross(repmat([B0(1) B0(2) B0(3)], size(normals_i,1),1),normals_i,2);
    B0xn_coeffs = proj(B0xn,Gamma_i,RWG_i);


    %% Non SS Computations 
    % Cross bilinear forms
    % grady G(x,y) . Vel(y)
    gradyGdotVely = @(x,y) 1/4/pi * (x-y)./vecnorm(x-y,2,2).^3 * (Vel(y))';
    gradxGdotVelx = @(x,y) -1/4/pi ./vecnorm(x-y,2,2).^3 .* dot(x-y,Vel(x),2);
    Gxy = @(x,y) 1/4/pi ./vecnorm(x-y,2,2);

    % ei partial derivative computation
    Aei1mat = integral(Gamma_i,Gamma_e,RWG_i,gradxGdotVelx,RWG_e);
    Aei1 = Psi_i' * Aei1mat * Psi_e;

    [Y_e,WY_e] = Gamma_e.qud;
    [X_i,WX_i] = Gamma_i.qud;
    NX = size(X_i,1);
    NY = size(Y_e,1);
    XX = repmat(X_i,NY,1); WWX = repmat(WX_i,NY,1); 
    YY = repelem(Y_e,NX,1); WWY = repelem(WY_e,NX,1);
    W = WWX .* WWY;
    VelXX = Vel(XX);
    DVel1XX = DVel{1}(XX);
    DVel2XX = DVel{2}(XX);
    DVel3XX = DVel{3}(XX);
    
    Psi_iXX = repmat(Psi_i_vals,NY,1);
    Psi_eYY = repelem(Psi_e_vals,NX,1);
    DVelPsi_iXX = [dot(DVel1XX,Psi_iXX,2) dot(DVel2XX,Psi_iXX,2) dot(DVel3XX,Psi_iXX,2)];
    Aei2 = 1/4/pi * sum(W.*(dot(Psi_eYY,DVelPsi_iXX,2)./(vecnorm(XX-YY,2,2))),1);

    % Cie computation
    % 1st term
    [X_e,WX_e] = Gamma_e.qud;
    [Y_i,WY_i] = Gamma_i.qud;
    NX = size(X_e,1);
    NY = size(Y_i,1);
    XX = repmat(X_e,NY,1); WWX = repmat(WX_e,NY,1); 
    YY = repelem(Y_i,NX,1); WWY = repelem(WY_i,NX,1);
    W = WWX .* WWY;

    nxgvals_iYY = repelem(nxgvals_i,NX,1);
    Psi_eXX = repmat(Psi_e_vals,NY,1);
    Ciekernel1 = 3/4/pi * (XX-YY).*dot(XX-YY,Vel(XX)-Vel(YY),2)./vecnorm(XX-YY,2,2).^5 ...
                -1/4/pi * (Vel(XX)-Vel(YY))./vecnorm(XX-YY,2,2).^3;

    Cie1 = sum(W.*dot(Psi_eXX,cross(Ciekernel1,nxgvals_iYY,2),2) ,1);

    % 2nd term
    DVel1YY = DVel{1}(YY);
    DVel2YY = DVel{2}(YY);
    DVel3YY = DVel{3}(YY);

    DVelnxgvals_iYY = [dot(DVel1YY,nxgvals_iYY,2) dot(DVel2YY,nxgvals_iYY,2) dot(DVel3YY,nxgvals_iYY,2)];
    % Gradx G
    Ciekernel2 = 1/4/pi * (YY-XX)./vecnorm(YY-XX,2,2).^3;
    Cie2 = sum(W.*dot(Psi_eXX,cross(Ciekernel2,DVelnxgvals_iYY,2) ,2) ,1);


    %% SS computations

    kernelA1 = @(x,y,z) dot(z,Vel(x) - Vel(y), 2)./(vecnorm(z,2,2).^3)/ (4*pi);

    kernelA2 = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
    % z := y-x, kernel is gradx G(x,y)
    kernelC1 = @(x,y,z) 1/(4*pi) * z./vecnorm(z,2,2).^3;
    
    % kernel for d/ds grad_xG(Ts(xh),Ts(yh))|_{s=0}
    kernelC3 = @(x,y,z) -3/(4*pi) * z .* dot(z,Vel(y)-Vel(x),2)./vecnorm(z,2,2).^5 + 1/(4*pi)*(Vel(y)-Vel(x))./vecnorm(z,2,2).^3;

    kernelN = kernelA1;


    Nelt_i = bndmesh_i.nelt;
    [ii,jj] = meshgrid(1:Nelt_i,1:Nelt_i);
    
    euler = parcluster('local');
    euler.NumWorkers = 5;
    saveProfile(euler);

    pool = euler.parpool(5);

    spmd
        if spmdIndex==1
            % partial derivative of b_A 
            A1mat_ii = panel_assembly(bndmesh_i,kernelA1,RWG_i,RWG_i,ii(:),jj(:));
            A1_ii = Psi_i' * A1mat_ii * Psi_i;

        elseif spmdIndex==2
            %g
            DVelRWG = RWG_i;
            DVelRWG.opr = 'Dvel[psi]';
            A2mat_ii = panel_assembly_shape_derivative(bndmesh_i,kernelA2,DVelRWG,RWG_i,ii(:),jj(:),Vel,DVel);
            A2_ii = 2* Psi_i' * A2mat_ii * Psi_i;

        elseif spmdIndex==3
            DVelRWG = RWG_i;
            DVelRWG.opr = 'Dvel[psi]';
            % Partial derivative of b_C
            C1mat_ii = panel_assembly_shape_derivative(bndmesh_i,kernelC1,DVelRWG,RWG_i,ii(:),jj(:),Vel,DVel);
            C1_ii = Psi_i' * C1mat_ii * g_i;
            C2_ii = g_i' * C1mat_ii * Psi_i ;% = C1?

        elseif spmdIndex==4
            % C3 (Is this way of evaluation okay?), z:= y-x
            C3mat_ii = panel_assembly(bndmesh_i,kernelC3,RWG_i,RWG_i,ii(:),jj(:));
            C3_ii = Psi_i' * C3mat_ii * g_i;

        elseif spmdIndex==5
            % Partial derivative of b_N
            Nmat_ii = panel_assembly_shape_derivative(bndmesh_i,kernelN,RWG_i.div,RWG_i.div,ii(:),jj(:),Vel,DVel);
            N_ii = -g_i' * Nmat_ii * g_i;

        end
    end
    
    % SS Based linear forms
    l1 = mu * jumpMuInv * Psi_i' * A1mat_ii{1} * B0xn_coeffs;

    NONEBEMSpace = RWG_i;
    NONEBEMSpace.opr = 'NONE';
    NONEBEMSpace.dir = B0;
    l2vec = panel_assembly_shape_derivative(bndmesh_i,kernelA2,NONEBEMSpace,RWG_i,ii(:),jj(:),Vel,DVel);
    l2 = mu * jumpMuInv * dot(l2vec,Psi_i);

    l3 = mu * jumpMuInv * B0xn_coeffs' * A2mat_ii{2} * Psi_i;

    l4 = mu0 * jumpMuInv * B0xn_coeffs' * C3mat_ii{4} * g_i;

    l5 = mu0 * jumpMuInv * B0xn_coeffs' * C1mat_ii{3} * g_i;

    l6vec = panel_assembly_shape_derivative(bndmesh_i,kernelC1,NONEBEMSpace,RWG_i,ii(:),jj(:),Vel,DVel);
    l6 = mu0 * jumpMuInv * dot(l6vec,g_i);

    [X_i,W_i] = Gamma_i.qud;
    DVelnxgvals_i = [dot(DVel{1}(X_i),nxgvals_i,2)  dot(DVel{2}(X_i),nxgvals_i,2) dot(DVel{3}(X_i),nxgvals_i,2)];
    l7integral = sum(W_i.* DVelnxgvals_i,1);
    l7 = -mu0 * jumpMuInv/2 * dot(B0,l7integral);

    r1 = -mu * jumpMuInv^2 /2 * B0xn_coeffs' * A1mat_ii{1} * B0xn_coeffs;
    r2 = -mu * jumpMuInv^2 /2 * dot(l2vec,B0xn_coeffs);
    r3 = r2;
    
    Vel_i = Vel(X_i);
    Veldotn_i = dot(Vel_i,normals_i,2);
    r4integral = sum(W_i.*Veldotn_i,1);
    r4 = -jumpMuInv/2 * norm(B0,2)^2 * r4integral;
    
    sd = -1/(2*mu0) * ( (1+mu/mu0) * (A1_ii{1}+A2_ii{2})...
                        +4 * (C1_ii{3}+C2_ii{3}+C3_ii{4})...
                        +(1+mu0/mu) * (N_ii{5}) ...
                        + 2 * (Aei1+Aei2)...
                        + 2 * (Cie1+Cie2))...
         +1/mu0 * (l1+l2+l3+l4+l5+l6+l7)...
         + r1 + r2 + r3 + r4;

    pool.delete();
end