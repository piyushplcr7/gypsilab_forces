function sd = SdBEMPMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel,DVel,mu0,H0,M)
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    % BEM Spaces
    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');
    P1_i = fem(bndmesh_i,'P1');

    g_i_vals = reconstruct(g_i,Gamma_i,P1_i);
    psi_e_vals = reconstruct(psi_e,Gamma_e,P0_e);

    normals_i = Gamma_i.qudNrm;

    %% Non SS Computations
    % Evaluating the Jacobian of the velocity field (row-wise) at gamma_i qud pts
    [X_i,W_i] = Gamma_i.qud;
    Mvals = M(X_i);
    Mdotn = dot(Mvals,normals_i,2);
    Mdotncoeffs = proj(Mdotn,Gamma_i,P0_i);

    Veli = Vel(X_i);
    DVel1i = DVel{1}(X_i);
    DVel2i = DVel{2}(X_i);
    DVel3i = DVel{3}(X_i);
    divVeli = DVel1i(:,1) + DVel2i(:,2) + DVel3i(:,3);
    divVelgi = divVeli.*g_i_vals;
    divVelgi_coeffs = proj(divVelgi,Gamma_i,P1_i);

    Kmatii = double_layer_laplace(Gamma_i,P0_i,P1_i);

    % Cross matrices
    dsVeiKernel = @(X,Y) 1/4/pi.* dot(Y-X,Vel(X),2) ./vecnorm(Y-X,2,2).^3;
    dsVei = integral(Gamma_i,Gamma_e,P0_i,dsVeiKernel,P0_e);

    dsKie1kernel = cell(3,1);
    dsKie1kernel{1} = @(X,Y) 3/4/pi * (X(:,1)-Y(1)) .* ((X-Y) * Vel(Y)')./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,Y,1) ./vecnorm(X-Y,2,2).^3;
    dsKie1kernel{2} = @(X,Y) 3/4/pi * (X(:,2)-Y(2)) .* ((X-Y) * Vel(Y)')./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,Y,2) ./vecnorm(X-Y,2,2).^3;
    dsKie1kernel{3} = @(X,Y) 3/4/pi * (X(:,3)-Y(3)) .* ((X-Y) * Vel(Y)')./vecnorm(X-Y,2,2).^5 - 1/4/pi * getFirstElem(Vel,Y,3) ./vecnorm(X-Y,2,2).^3;
    
    dsKie1 = integral(Gamma_e,Gamma_i,P0_e,dsKie1kernel,ntimes(P1_i));

    % dsKie2 computation explicit
    [Y,WY] = Gamma_i.qud;
    [X,WX] = Gamma_e.qud;

    NX = size(X,1); NY = size(Y,1);

    XX = repelem(X,NY,1); YY = repmat(Y,NX,1);
    W = repelem(WX,NY,1).*repmat(WY,NX,1);

    g_i_vals_YY = repmat(g_i_vals,NX,1);
    psi_e_XX = repelem(psi_e_vals,NY,1);
    normalsi_YY = repmat(normals_i,NX,1);
    DVel1iYY = repmat(DVel1i,NX,1);
    DVel2iYY = repmat(DVel2i,NX,1);
    DVel3iYY = repmat(DVel3i,NX,1);
    divVeliYY = repmat(divVeli,NX,1);

    gradyGXXYY = 1/4/pi * (XX-YY)./vecnorm(XX-YY,2,2).^3;
    DVelgradyGXXYY = [dot(DVel1iYY,gradyGXXYY,2) dot(DVel2iYY,gradyGXXYY,2) dot(DVel3iYY,gradyGXXYY,2)];

    dskie2kernel = divVeliYY.*dot(gradyGXXYY,normalsi_YY,2) - dot(normalsi_YY,DVelgradyGXXYY,2);
    dsKie2 = sum(W.*dskie2kernel.*g_i_vals_YY.*psi_e_XX,1);


    %% SS Computations
    Nelt_i = bndmesh_i.nelt;

    [ii,jj] = meshgrid(1:Nelt_i,1:Nelt_i);

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
            kerneloldmat_P0_P0_ii = panel_assembly(bndmesh_i,kernelold,P0_i,P0_i,ii(:),jj(:));
            % Partial derivative of bv(psi,psi)
            dbv_dsii = psi_i' * kerneloldmat_P0_P0_ii * psi_i;

        elseif spmdIndex==2
            kerneloldmat_nxgradP1_nxgradP1ii = panel_assembly(bndmesh_i,kernelold,P1_i.nxgrad,P1_i.nxgrad,ii(:),jj(:));

        elseif spmdIndex==3
            kernelintegrablematii = panel_assembly(bndmesh_i,kernelintegrable,ntimes(P1_i),P0_i,ii(:),jj(:));

        elseif spmdIndex==4
            % Combination kernel that cancels singularity
            combkernelmatii = panel_assembly(bndmesh_i,combkernel,ntimes(P1_i),P0_i,ii(:),jj(:));

        elseif spmdIndex==5
             % Partial derivative of bw(g,g)
            SL_Dvelnxgrad_nxgradii = panel_assembly_shape_derivative(bndmesh_i,KV,P1_i.nxgrad,P1_i.nxgrad,ii(:),jj(:),Vel,DVel);
        end
    end

    dbk_dsii = psi_i' * Kmatii * divVelgi_coeffs + psi_i' * (kernelintegrablematii{3} -combkernelmatii{4}) * g_i;
    dbw_dsii = g_i' * ( kerneloldmat_nxgradP1_nxgradP1ii{2} + 2 * SL_Dvelnxgrad_nxgradii{5}) * g_i;

    % Linear form
    l1 = -psi_i' * kerneloldmat_P0_P0_ii{1} * Mdotncoeffs;
    l23 = -(Mdotncoeffs' * Kmatii * divVelgi_coeffs + Mdotncoeffs' * (kernelintegrablematii{3} -combkernelmatii{4}) * g_i);

    % Remaining terms
    % DVelMi = [dot(DVel1i,Mvals,2) dot(DVel2i,Mvals,2) dot(DVel3i,Mvals,2)];
    % DVelMidotH0 = DVelMi * H0';
    % r1 = mu0 * sum(W_i.*DVelMidotH0,1);
    r1 = mu0 * sum(W_i.*Mdotn.*(Veli * H0'),1);
    r2 = -mu0/2 * Mdotncoeffs' * kerneloldmat_P0_P0_ii{1} * Mdotncoeffs;


    sd = -mu0/2*( 2 * dbv_dsii{1} + 4 * dbk_dsii -2 * dbw_dsii...
                  + 2 * psi_i' * dsVei * psi_e + 2 * psi_e' * dsKie1 * g_i + 2 * dsKie2)...
         +mu0 * (l1+l23)...
         + r1 + r2;

    pool.delete();
end