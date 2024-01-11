function sd = sdBEMTpDielSC(bndmesh,epsilon,epsilon0,Tnu,Tdu,rho,omega_src,Vel,DVel)
    P0 = fem(bndmesh,'P0');
    P1 = fem(bndmesh,'P1');

    nxgradP1 = P1.nxgrad;

    Gamma = dom(bndmesh,3);
    [Xgypsi,~] = Gamma.qud;
    normals = Gamma.qudNrm;

    % Reconstructing the trace of u- > g
    g = reconstruct(Tdu,Gamma,P1);

    %% Non SS computations
    
    % Evaluating the velocity field at quadrature points
    % Vels = Vel(Xgypsi);
    % Evaluating the Jacobian of the velocity field (row-wise) at qud pts
    DVel1 = DVel{1}(Xgypsi);
    DVel2 = DVel{2}(Xgypsi);
    DVel3 = DVel{3}(Xgypsi);
    divVel = DVel1(:,1) + DVel2(:,2) + DVel3(:,3);

    % partial derivative of bk(g,psi)
    divVelg = divVel.*g;
    divVelg_coeffs = proj(divVelg,Gamma,P1);
    Kmat = double_layer_laplace(Gamma,P0,P1);

    %% SS Computations

    Nelt = bndmesh.nelt;

    [ii,jj] = meshgrid(1:Nelt,1:Nelt);

    % Kernel gradxG.vel(x) + gradyG.vel(y), z:= y-x
    kernelold = @(x,y,z) dot(z,Vel(x) - Vel(y), 2)./(vecnorm(z,2,2).^3)/ (4*pi);

    kernelintegrable = @(x,y,z) 3/(4*pi)* dot(z,Vel(y) - Vel(x),2) .*z ./vecnorm(z,2,2).^5;

    combkernel = @(x,y,z) 1/(4*pi) * ( -[ dot(DVel{1}(y),z,2) dot(DVel{2}(y),z,2) dot(DVel{3}(y),z,2) ] + Vel(y) - Vel(x) )./vecnorm(z,2,2).^3;
    
    KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;

    euler = parcluster('local');
    euler.NumWorkers = 5;
    saveProfile(euler);

    pool = euler.parpool(5);

    spmd
        if spmdIndex==1
            kerneloldmat_P0_P0 = panel_assembly(bndmesh,kernelold,P0,P0,ii(:),jj(:));
            % Partial derivative of bv(psi,psi)
            dbv_ds = Tnu' * kerneloldmat_P0_P0 * Tnu;

        elseif spmdIndex==2
            kerneloldmat_nxgradP1_nxgradP1 = panel_assembly(bndmesh,kernelold,nxgradP1,nxgradP1,ii(:),jj(:));

        elseif spmdIndex==3
            kernelintegrablemat = panel_assembly(bndmesh,kernelintegrable,ntimes(P1),P0,ii(:),jj(:));

        elseif spmdIndex==4
            % Combination kernel that cancels singularity
            combkernelmat = panel_assembly(bndmesh,combkernel,ntimes(P1),P0,ii(:),jj(:));

        elseif spmdIndex==5
             % Partial derivative of bw(g,g)
            SL_Dvelnxgrad_nxgrad = panel_assembly_shape_derivative(bndmesh,KV,nxgradP1,nxgradP1,ii(:),jj(:),Vel,DVel);
        end
    end
    
    dbk_ds_reduced = Tnu' * (kernelintegrablemat{3} -combkernelmat{4}) * Tdu;
    dbk_ds = Tnu' * Kmat * divVelg_coeffs + Tnu' * (kernelintegrablemat{3} -combkernelmat{4}) * Tdu;
    dbw_ds = Tdu' * ( kerneloldmat_nxgradP1_nxgradP1{2} + 2 * SL_Dvelnxgrad_nxgrad{5}) * Tdu;

    %% Copied from Const Vel
    psi_vals = reconstruct(Tnu,Gamma,P0);
    g_vals = reconstruct(Tdu,Gamma,P1);

    [Y,WY] = omega_src.qud;
    [X,WX] = Gamma.qud;
    NY = size(Y,1);
    NX = size(X,1);

    YY = repmat(Y,NX,1);
    XX = repelem(X,NY,1);

    W = repmat(WY,NX,1).*repelem(WX,NY,1);

    gradxGdotvelx = 1/4/pi * dot(YY-XX,Vel(XX),2)./vecnorm(XX-YY,2,2).^3;
    rhoYY = rho(YY);
    psiXX = repelem(psi_vals,NY,1);
    gXX = repelem(g_vals,NY,1);

    l1 = - sum(W.* gradxGdotvelx.*rhoYY.*psiXX ,1);

    normals_XX = repelem(normals,NY,1);

    t21kernel = 3/4/pi * dot(XX-YY,normals_XX,2).*dot(XX-YY,Vel(XX),2)./vecnorm(XX-YY,2,2).^5 - 1/4/pi * dot(Vel(XX),normals_XX,2)./vecnorm(XX-YY,2,2).^3;

    divVeln = divVel.*normals;
    
    DVel1XX = DVel{1}(XX);
    DVel2XX = DVel{2}(XX);
    DVel3XX = DVel{3}(XX);

    gradxGXXYY = 1/4/pi * (YY-XX)./vecnorm(XX-YY,2,2).^3;

    DVelgradxGXX = [dot(DVel1XX,gradxGXXYY,2) dot(DVel2XX,gradxGXXYY,2) dot(DVel3XX,gradxGXXYY,2)];

    divVelnXX = repelem(divVeln,NY,1);

    l3 = sum(W.*(dot(gradxGXXYY,divVelnXX,2) - dot(normals_XX,DVelgradxGXX,2)).*rhoYY.*gXX,1);


    l2 = sum(W.*t21kernel.*rhoYY.*gXX,1);


    %% Final SD
    sd = epsilon/2 * ( (1+epsilon0/epsilon) * dbv_ds{1}...
                       -4 * dbk_ds...
                       -(1+epsilon/epsilon0) * dbw_ds)...
                      +l1 + l2 + l3;

    pool.delete();
end