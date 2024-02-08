function sd = SdBEMLMCFSP_dualnorm_ALT(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,mu0,mu,H0,abc_alpha,order)
    Nfields = size(abc_alpha,1);
    % order = 3;
    Gamma_i = dom(bndmesh_i,order);
    Gamma_e = dom(bndmesh_e,order);

    % BEM Spaces
    P0_i = fem(bndmesh_i,'P0');
    P1_i = fem(bndmesh_i,'P1');
    P1_e = fem(bndmesh_e,'P1');

    g_i_vals = reconstruct(g_i,Gamma_i,P1_i);

    % Dirichlet trace at the outer boundary
    [X_e,~] = Gamma_e.qud;
    g_e_vals = X_e * H0';
    g_e = proj(g_e_vals,Gamma_e,P1_e);


    %% Dispatching SS computations to GPU 
    Tnu = psi_i;
    Tdu = g_i;
    bndmesh = bndmesh_i;
    P1 = P1_i;
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
    
    ptxFilePath = 'SdBEMLMCFSP_dualnorm_GPU_ALT.ptx';
    cuFilePath = 'SdBEMLMCFSP_dualnorm_GPU_ALT.cu';

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
    H0dotn_coeffs = 0;
    H0dotn_coeffs_gpu = gpuArray(full(H0dotn_coeffs));
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
    Tdu_gpu, Tnu_gpu,H0dotn_coeffs_gpu,...
    Elements_gpu,Vertices_gpu,Normals_gpu,Areas_gpu,...
    elt2dof_gpu,elt2dof_gpu,...
    TrialSpace_gpu,TestSpace_gpu,TrialOperator_gpu,TestOperator_gpu,...
    size(P1.rsf,1),size(P1.rsf,1),...
    abc_alpha_gpu,Nfields,...
    permI_gpu,permJ_gpu); 

    %% Non SS Computations
    [X_i,W_i] = Gamma_i.qud;
    Kmatii = double_layer_laplace(Gamma_i,P0_i,P1_i);
    sdnongpu = zeros(Nfields,1);
    dbk_dsii = zeros(Nfields,1);
    dsVie = zeros(Nfields,1);
    dsKie1 = zeros(Nfields,1);
    dsKie2 = zeros(Nfields,1);
    dsl1 = zeros(Nfields,1);
    dsl2 = zeros(Nfields,1);
    dsl3 = zeros(Nfields,1);

    for fieldID = 1:Nfields
        a = abc_alpha(fieldID,1);
        b = abc_alpha(fieldID,2);
        c = abc_alpha(fieldID,3);
        alpha = abc_alpha(fieldID,4);
        [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);
        DVel1i = DVel{1}(X_i);
        DVel2i = DVel{2}(X_i);
        DVel3i = DVel{3}(X_i);
        divVeli = DVel1i(:,1) + DVel2i(:,2) + DVel3i(:,3);
        divVelgi = divVeli.*g_i_vals;
        divVelgi_coeffs = proj(divVelgi,Gamma_i,P1_i);

        dbk_dsii(fieldID) = psi_i' * Kmatii * divVelgi_coeffs;

        dsVie(fieldID) = dsVieFun(bndmesh_i,bndmesh_e,psi_i,psi_e,Vel,DVel);
        dsKie1(fieldID) = dsKieFun1(bndmesh_i,bndmesh_e,g_i,psi_e,Vel,DVel);
        dsKie2(fieldID) = dsKieFun2(bndmesh_i,bndmesh_e,g_i,psi_e,Vel,DVel);

        dsl1(fieldID) = dsl1_ALT(bndmesh_i,bndmesh_e,g_e,psi_i,Vel,DVel);
        dsl2(fieldID) = dsl2_ALT(bndmesh_i,bndmesh_e,g_e,g_i,Vel,DVel);
        dsl3(fieldID) = dsl3_ALT(bndmesh_i,bndmesh_e,g_e,g_i,Vel,DVel);
        
    end

    sdnongpu = -mu0/2*( 4 * dbk_dsii ...
                  + 2 * dsVie + 2 * dsKie1 + 2 * dsKie2)...
         +mu0 * (dsl1+dsl2+dsl3);

    sd = sdnongpu + shapeDerivative_gpu;

end