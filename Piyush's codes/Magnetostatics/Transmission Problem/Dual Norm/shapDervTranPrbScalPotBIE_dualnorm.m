% Full shape Derivative

function sd_gpu_full = shapDervTranPrbScalPotBIE_dualnorm(bndmesh,Tdu,Tnu,J,omega_src,mu0,mu,abc_alpha)
    Nfields = size(abc_alpha,1);
    % BEM Spaces
    P0 = fem(bndmesh,'P0');
    P1 = fem(bndmesh,'P1');
    gradP1 = P1.grad;
    nxgradP1 = P1.nxgrad;

    Gamma = dom(bndmesh,3);

    [Xgypsi,Wgypsi] = Gamma.qud;
    normals = Gamma.qudNrm;

    % Reconstructing the Neumann trace psi
    psi = reconstruct(Tnu,Gamma,P0);
    % Reconstructing the surface gradient of u
    gradTu = reconstruct(Tdu,Gamma,gradP1);
    % Reconstructing the trace of u- > g
    g = reconstruct(Tdu,Gamma,P1);

    HJ = compute_vecpot_curl(J,omega_src,Xgypsi);
    DHJ = compute_vecpot_D_curl(J,omega_src,Xgypsi);
    HJn = dot(normals,HJ,2);
    HJn_coeffs = proj(HJn,Gamma,P0);

    jumpMu = mu0-mu;
    %% Full shape derivative

    %% Non SS computations

    % partial derivative of bk(g,psi)
    Kmat = double_layer_laplace(Gamma,P0,P1);
    slmat = single_layer(Gamma,P0,P0);

    %% Dispatching SS computations to GPU 
    [~,elt2dof] = P1.dof;
    % Creating the interaction matrix
    elts = 1:bndmesh.nelt;
    elts = elts';
    elts = repelem(elts,3);
    vtcs = reshape(bndmesh.elt',[bndmesh.nelt*3 1]);
    
    Eltmat = sparse(elts,vtcs,ones(size(vtcs)),bndmesh.nelt,bndmesh.nvtx);
    
    Intmat = Eltmat * Eltmat';
    [I0,J0] = find(Intmat == 0);
    [I1,J1] = find(Intmat == 1);
    [I2,J2] = find(Intmat == 2);
    [I3,J3] = find(Intmat == 3);
    Ivec = [I0;I1;I2;I3]-1;
%     Ivec = cast(Ivec,'int32');
    Ivec = cast(Ivec,'uint16');
    Jvec = [J0;J1;J2;J3]-1;
%     Jvec = cast(Jvec,'int32');
    Jvec = cast(Jvec,'uint16');
    relation = [0*ones(size(I0,1),1);...
                1*ones(size(I1,1),1);...
                2*ones(size(I2,1),1);...
                3*ones(size(I3,1),1)];
    relation = cast(relation,'int32');

    % Creating vectors permI and permJ

    % For non-interacting elements, perm is trivial
    permI0 = repmat([0 1 2],size(I0,1),1);
    permJ0 = permI0;

    % Elements with interaction == 1
    eltI1 = bndmesh.elt(I1,:);
    eltJ1 = bndmesh.elt(J1,:);
    [intersection1,diffI1,diffJ1] = rowWiseIntersectionDiff(eltI1,eltJ1);
    assert(size(intersection1,2) == 1);
    assert(size(diffI1,2) == 2);
    assert(size(diffJ1,2) == 2);
    permI1 = findPermVectorized([intersection1 diffI1],eltI1);
    permJ1 = findPermVectorized([intersection1 diffJ1],eltJ1);

    % Elements with interaction == 2
    eltI2 = bndmesh.elt(I2,:);
    eltJ2 = bndmesh.elt(J2,:);
    [intersection2,diffI2,diffJ2] = rowWiseIntersectionDiff(eltI2,eltJ2);
    assert(size(intersection2,2) == 2);
    assert(size(diffI2,2) == 1);
    assert(size(diffJ2,2) == 1);
    permI2 = findPermVectorized([intersection2 diffI2],eltI2);
    permJ2 = findPermVectorized([intersection2 diffJ2],eltJ2);

    % Elements with interaction == 3, perm is trivial
    permI3 = repmat([0 1 2],size(I3,1),1);
    permJ3 = permI3;

    permI = [permI0;(permI1-1);(permI2-1);permI3];
    permJ = [permJ0;(permJ1-1);(permJ2-1);permJ3];

    permI_gpu = gpuArray(cast(permI','int32'));
    permJ_gpu = gpuArray(cast(permJ','int32'));

    % PTX and CUDA files
    ptxFilePath = 'shapDervTranPrbScalPotBIE_dualnorm_GPU.ptx';
    cuFilePath = 'shapDervTranPrbScalPotBIE_dualnorm_GPU.cu';

    % Quadrature to be passed to the GPU
    order = 5;
    [X, W] = quad4D(order);

    % Kernel Object
    kernel = parallel.gpu.CUDAKernel(ptxFilePath, cuFilePath);
    gridDim = [300, 1, 1];
    blockDim = [32, 1, 1]; % 32 is the SIMD Width or wrap size
    Nthreads = prod(blockDim) * prod(gridDim);
    kernel.GridSize = gridDim;
    kernel.ThreadBlockSize = blockDim; 

    % Transferring Interactions and Relations to GPU
    Ivec_gpu = gpuArray(Ivec);
    Jvec_gpu = gpuArray(Jvec);
    relation_gpu = gpuArray(relation);

    % Transferring Quadrature to GPU
    W0_gpu = gpuArray(W{4});
    W1_gpu = gpuArray(W{3});
    W2_gpu = gpuArray(W{2});
    W3_gpu = gpuArray(W{1});
    X0 = X{4}; X0(:,1) = X0(:,1)-X0(:,2); X0(:,3) = X0(:,3)-X0(:,4);
    X0_gpu = gpuArray(X0');
    X1 = X{3}; X1(:,1) = X1(:,1)-X1(:,2); X1(:,3) = X1(:,3)-X1(:,4);
    X1_gpu = gpuArray(X1');
    X2 = X{2}; X2(:,1) = X2(:,1)-X2(:,2); X2(:,3) = X2(:,3)-X2(:,4);
    X2_gpu = gpuArray(X2');
    X3 = X{1}; X3(:,1) = X3(:,1)-X3(:,2); X3(:,3) = X3(:,3)-X3(:,4);
    X3_gpu = gpuArray(X3');

    % Output variables GPU
    dbv_ds_gpu = gpuArray.zeros(1,1);
    dbw_ds_gpu = gpuArray.zeros(1,1);
    dbk_ds_reduced_gpu = gpuArray.zeros(1,1); % excludes the Kmat part which is computed by the CPU
    shapeDerivative_gpu = gpuArray.zeros(Nfields,1);
    
    % More input variables to GPU
    Tdu_gpu = gpuArray(Tdu);
    Tnu_gpu = gpuArray(Tnu);
    HJn_coeffs_gpu = gpuArray(full(HJn_coeffs));
    Elements = cast(bndmesh.elt-1,'int32');
    Elements_gpu = gpuArray(Elements');
    Vertices_gpu = gpuArray(bndmesh.vtx');
    Normals_gpu = gpuArray(bndmesh.nrm');
    Areas_gpu = gpuArray(bndmesh.ndv);
    zeroIdxelt2dof = cast(elt2dof-1,'int32');
    elt2dof_gpu = gpuArray(zeroIdxelt2dof');
    TrialSpace_gpu = 0;
    TestSpace_gpu = 0;
    TrialOperator_gpu = 0;
    TestOperator_gpu = 0;
    abc_alpha_gpu = gpuArray(cast(abc_alpha','int32'));

    % Launching kernel
    [dbv_ds_gpu,dbw_ds_gpu,dbk_ds_reduced_gpu, shapeDerivative_gpu]= feval(kernel,...
    0,0,bndmesh.nelt,bndmesh.nvtx,bndmesh.nelt^2,...
    Nthreads,Ivec_gpu,Jvec_gpu,relation_gpu,...
    W0_gpu,X0_gpu,size(X0,1),...
    W1_gpu,X1_gpu,size(X1,1),...
    W2_gpu,X2_gpu,size(X2,1),...
    W3_gpu,X3_gpu,size(X3,1),...
    dbv_ds_gpu,dbw_ds_gpu,dbk_ds_reduced_gpu,...
    shapeDerivative_gpu,...
    mu0,mu,...
    Tdu_gpu, Tnu_gpu, HJn_coeffs_gpu,...
    Elements_gpu,Vertices_gpu,Normals_gpu,Areas_gpu,...
    elt2dof_gpu,elt2dof_gpu,...
    TrialSpace_gpu,TestSpace_gpu,TrialOperator_gpu,TestOperator_gpu,...
    size(P1.rsf,1),size(P1.rsf,1),...
    abc_alpha_gpu,Nfields,...
    permI_gpu,permJ_gpu); 


    %% SS Computations

%     Nelt = bndmesh.nelt;

%     [ii,jj] = meshgrid(1:Nelt,1:Nelt);
% 
%     % Kernel gradxG.vel(x) + gradyG.vel(y), z:= y-x
%     kernelold = @(x,y,z) dot(z,Vel(x) - Vel(y), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
% 
%     kernelintegrable = @(x,y,z) 3/(4*pi)* dot(z,Vel(y) - Vel(x),2) .*z ./vecnorm(z,2,2).^5;
% 
%     combkernel = @(x,y,z) 1/(4*pi) * ( -[ dot(DVel{1}(y),z,2) dot(DVel{2}(y),z,2) dot(DVel{3}(y),z,2) ] + Vel(y) - Vel(x) )./vecnorm(z,2,2).^3;
%     
%     KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
% 
%     euler = parcluster('local');
%     euler.NumWorkers = 5;
%     saveProfile(euler);
% 
%     pool = euler.parpool(5);

%     spmd
%         if spmdIndex==1
%             kerneloldmat_P0_P0 = panel_assembly(bndmesh,kernelold,P0,P0,ii(:),jj(:));
%             % Partial derivative of bv(psi,psi)
%             dbv_ds = Tnu' * kerneloldmat_P0_P0 * Tnu;
% 
%         elseif spmdIndex==2
%             kerneloldmat_nxgradP1_nxgradP1 = panel_assembly(bndmesh,kernelold,nxgradP1,nxgradP1,ii(:),jj(:));
% 
%         elseif spmdIndex==3
%             kernelintegrablemat = panel_assembly(bndmesh,kernelintegrable,ntimes(P1),P0,ii(:),jj(:));
% 
%         elseif spmdIndex==4
%             % Combination kernel that cancels singularity
%             combkernelmat = panel_assembly(bndmesh,combkernel,ntimes(P1),P0,ii(:),jj(:));
% 
%         elseif spmdIndex==5
%              % Partial derivative of bw(g,g)
%             SL_Dvelnxgrad_nxgrad = panel_assembly_shape_derivative(bndmesh,KV,nxgradP1,nxgradP1,ii(:),jj(:),Vel,DVel);
%         end
%     end

%     dbk_ds = Tnu' * Kmat * divVelg_coeffs + Tnu' * (kernelintegrablemat{3} -combkernelmat{4}) * Tdu;
%     dbk_ds_reduced = Tnu' * (kernelintegrablemat{3} -combkernelmat{4}) * Tdu;
% 
%    
%     dbw_ds = Tdu' * ( kerneloldmat_nxgradP1_nxgradP1{2} + 2 * SL_Dvelnxgrad_nxgrad{5}) * Tdu;
    
    
%     dl2_ds = zeros(Nfields,1);
%     remterm1 = zeros(Nfields,1);
    elseterms = zeros(Nfields,1);

    for fieldID = 1:Nfields
        a = abc_alpha(fieldID,1);
        b = abc_alpha(fieldID,2);
        c = abc_alpha(fieldID,3);
        alpha = abc_alpha(fieldID,4);
        [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);

        % Evaluating the velocity field at quadrature points
        Vels = Vel(Xgypsi);
        % Evaluating the Jacobian of the velocity field (row-wise) at qud pts
        DVel1 = DVel{1}(Xgypsi);
        DVel2 = DVel{2}(Xgypsi);
        DVel3 = DVel{3}(Xgypsi);
        divVel = DVel1(:,1) + DVel2(:,2) + DVel3(:,3);
    
        divVelg = divVel.*g;
        divVelg_coeffs = proj(divVelg,Gamma,P1);
    
        % Projecting the complex integrand to P0 space
        DHJV = [dot(DHJ{1},Vels,2) dot(DHJ{2},Vels,2) dot(DHJ{3},Vels,2)];
        DVHJ = [dot(DVel1,HJ,2) dot(DVel2,HJ,2) dot(DVel3,HJ,2)];
        compl_integrand = dot(normals,(DHJV - DVHJ + HJ.*divVel),2);
        compl_integrand_coeffs = proj(compl_integrand,Gamma,P0);

        Vn = dot(Vels,normals,2);

        dl2_ds = sum(Wgypsi.*compl_integrand.*g,1);
        remterm1 = -0.5 * jumpMu * sum(Wgypsi.* dot(HJ,HJ,2) .* Vn,1);

        elseterms(fieldID) = mu0/2 * 4 * Tnu' * Kmat * divVelg_coeffs...
                 -mu0 * ( jumpMu/mu * Tnu' * slmat * compl_integrand_coeffs...
                 + jumpMu/2/mu0 * dl2_ds...
                 -jumpMu/mu0 * (HJn_coeffs' * Kmat * divVelg_coeffs...
                                + compl_integrand_coeffs' * Kmat * Tdu))...
                 + remterm1...
                 + -jumpMu^2/2/mu *HJn_coeffs' * 2* slmat * compl_integrand_coeffs;
    end
    
    
    
    

    % partial derivative of l1(psi)
%     dl1_ds = Tnu' * (kerneloldmat_P0_P0{1} * HJn_coeffs + slmat * compl_integrand_coeffs);

    % partial derivative of l2(g)
%     dl2_ds = sum(Wgypsi.*compl_integrand.*g,1);

    % partial derivative of l3(g)
%     dl3_ds = HJn_coeffs' * Kmat * divVelg_coeffs + HJn_coeffs' * (kernelintegrablemat{3} -combkernelmat{4}) * Tdu...
%             + compl_integrand_coeffs' * Kmat * Tdu;
    
    % Remaining terms
    
%     remterm1 = -0.5 * jumpMu * sum(Wgypsi.* dot(HJ,HJ,2) .* Vn,1);
%     remterms = -jumpMu^2/2/mu *HJn_coeffs' * (kerneloldmat_P0_P0{1} * HJn_coeffs + 2* slmat * compl_integrand_coeffs);

%     sd = mu0/2 * ( - (1+mu0/mu) * dbv_ds{1} ...
%                    + 4 * dbk_ds...
%                    + (1+mu/mu0) * dbw_ds)...
%          -mu0 * ( jumpMu/mu * dl1_ds ...
%                   + jumpMu/2/mu0 * dl2_ds ...
%                   -jumpMu/mu0 * dl3_ds)...
%                   +remterm1 + remterms;

%     sd_gpu = mu0/2 *( - (1+mu0/mu) * dbv_ds{1} ...
%                    + 4 * dbk_ds_reduced...
%                    + (1+mu/mu0) * dbw_ds)...
%                    -mu0 * (jumpMu/mu * Tnu' * kerneloldmat_P0_P0{1} * HJn_coeffs...
%                    -jumpMu/mu0 * HJn_coeffs' * (kernelintegrablemat{3} -combkernelmat{4}) * Tdu)...
%                    -jumpMu^2/2/mu *HJn_coeffs' * kerneloldmat_P0_P0{1} * HJn_coeffs;

%     sd_gpu_full = shapeDerivative_gpu...
%                  +mu0/2 * 4 * Tnu' * Kmat * divVelg_coeffs...
%                  -mu0 * ( jumpMu/mu * Tnu' * slmat * compl_integrand_coeffs...
%                  + jumpMu/2/mu0 * dl2_ds...
%                  -jumpMu/mu0 * (HJn_coeffs' * Kmat * divVelg_coeffs...
%                                 + compl_integrand_coeffs' * Kmat * Tdu))...
%                  + remterm1...
%                  + -jumpMu^2/2/mu *HJn_coeffs' * 2* slmat * compl_integrand_coeffs;

sd_gpu_full = shapeDerivative_gpu + elseterms;
    
    % sd = mu0/2 * ( - (1+mu0/mu) * dbv_ds ...
    %                + 4 * dbk_ds...
    %                + (1+mu/mu0) * dbw_ds)...
    %      -mu0 * ( jumpMu/mu * dl1_ds ...
    %               + jumpMu/2/mu0 * dl2_ds ...
    %               -jumpMu/mu0 * dl3_ds)...
    %               +remterm1 + remterms;

%     pool.delete();

end
