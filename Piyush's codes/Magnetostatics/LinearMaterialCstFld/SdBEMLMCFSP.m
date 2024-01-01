function sd = SdBEMLMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel,DVel,mu0,mu,H0)
    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);

    jumpMu = mu0 - mu;

    % BEM Spaces
    P0_i = fem(bndmesh_i,'P0');
    P0_e = fem(bndmesh_e,'P0');
    P1_i = fem(bndmesh_i,'P1');
    P1_e = fem(bndmesh_e,'P1');

    psi_i_vals = reconstruct(psi_i,Gamma_i,P0_i);
    g_i_vals = reconstruct(g_i,Gamma_i,P1_i);
    psi_e_vals = reconstruct(psi_e,Gamma_e,P0_e);

    %% Non SS Computations
    % Evaluating the Jacobian of the velocity field (row-wise) at gamma_i qud pts
    [X_i,W_i] = Gamma_i.qud;
    Veli = Vel(X_i);
    DVel1i = DVel{1}(X_i);
    DVel2i = DVel{2}(X_i);
    DVel3i = DVel{3}(X_i);
    divVeli = DVel1i(:,1) + DVel2i(:,2) + DVel3i(:,3);
    divVelgi = divVeli.*g_i_vals;
    divVelgi_coeffs = proj(divVelgi,Gamma_i,P1_i);

    Kmatii = double_layer_laplace(Gamma_i,P0_i,P1_i);
    Vii = single_layer(Gamma_i,P0_i,P0_i);

    % Cross matrices
    Vei = single_layer_cross(Gamma_i,Gamma_e,P0_i,P0_e);
    Kie = double_layer_laplace_cross(Gamma_e,Gamma_i,P0_e,P1_i);
    Vee = single_layer(Gamma_e,P0_e,P0_e);

    % H0.n coefficients
    normals_i = Gamma_i.qudNrm;
    H0extended = repmat(H0,size(normals_i,1),1);
    H0dotn_vals = dot(H0extended,normals_i,2);
    H0dotn_coeffs = proj(H0dotn_vals,Gamma_i,P0_i);

    % Complex integrand
    DVelH0 = [dot(DVel1i,H0extended,2) dot(DVel2i,H0extended,2) dot(DVel3i,H0extended,2)];
    compl_integrand = dot(normals_i,H0extended.*divVeli - DVelH0,2);
    compl_integrand_coeffs = proj(compl_integrand,Gamma_i,P0_i);

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
    l1 = jumpMu/mu * psi_i' * kerneloldmat_P0_P0_ii{1} * H0dotn_coeffs;
    l2 = jumpMu/mu * psi_i' * Vii * compl_integrand_coeffs;
    l3 = -jumpMu/2/mu0 * g_i' *  mass_matrix(Gamma_i,P1_i,P0_i) * compl_integrand_coeffs;
    l45 = jumpMu/mu0 * (H0dotn_coeffs' * Kmatii * divVelgi_coeffs + H0dotn_coeffs' * (kernelintegrablematii{3} -combkernelmatii{4}) * g_i);
    l6 = jumpMu/mu * compl_integrand_coeffs' * Kmatii * g_i;

    % Remaining terms
    r1 = -jumpMu^2/2/mu * H0dotn_coeffs' * kerneloldmat_P0_P0_ii{1} * H0dotn_coeffs;
    r2 = -jumpMu^2/2/mu * H0dotn_coeffs' * Vii * compl_integrand_coeffs;
    r3 = -jumpMu^2/2/mu * compl_integrand_coeffs' * Vii * H0dotn_coeffs;
    r4 = -jumpMu/2 * norm(H0)^2 * sum(W_i.*dot(Veli,normals_i,2),1);


    sd = -mu0/2*( (1+mu0/mu) * dbv_dsii{1} + 4 * dbk_dsii -(1+mu/mu0) * dbw_dsii...
                  + 2 * psi_i' * Vei * psi_e + 2 * psi_e' * Kie * g_i...
                  + psi_e' * Vee * psi_e)...
         +mu0 * (l1+l2+l3+l45+l6)...
         + r1 + r2 + r3 + r4;

    pool.delete();
end