function sd = SdBemLMCFVP_dualnorm(bndmesh_i,bndmesh_e,Psi_i,g_i,Psi_e,mu0,mu,B0,abc_alpha)
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

    Nfields = size(abc_alpha,1);

    %% Dispatching SS computations to GPU
    RWG = RWG_i;
    bndmesh = bndmesh_i;
    TnA = Psi_i;
    TdA = g_i;

    [~,elt2dof] = RWG.dof;
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
    ptxFilePath = 'SdBEMLMCFVP_dualnorm_GPU.ptx';
    cuFilePath = 'SdBEMLMCFVP_dualnorm_GPU.cu';

    % Quadrature to be passed to the GPU
    order = 5;
    [X, W] = quad4D(order);
    
    % Kernel Object
    kernel = parallel.gpu.CUDAKernel(ptxFilePath, cuFilePath);
    gridDim = [20000, 1, 1];
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

    % More variables for GPU
    A1_gpu = gpuArray.zeros(Nfields,1);
    A2_gpu = gpuArray.zeros(Nfields,1);
    C1_gpu = gpuArray.zeros(Nfields,1);
    C2_gpu = gpuArray.zeros(Nfields,1);
    C3_gpu = gpuArray.zeros(Nfields,1);
    N_gpu = gpuArray.zeros(Nfields,1);
    l1_gpu = gpuArray.zeros(Nfields,1);
    l2_gpu = gpuArray.zeros(Nfields,1);
    l3_gpu = gpuArray.zeros(Nfields,1);
    l4_gpu = gpuArray.zeros(Nfields,1);
    l5_gpu = gpuArray.zeros(Nfields,1);
    l6_gpu = gpuArray.zeros(Nfields,1);
    r1_gpu = gpuArray.zeros(Nfields,1);
    r2r3_gpu = gpuArray.zeros(Nfields,1);
    l2vec_gpu = gpuArray.zeros(size(TnA,1),1);

    shapeDerivative_gpu = gpuArray.zeros(Nfields,1);
    TnA_gpu = gpuArray(full(TnA));
    TdA_gpu = gpuArray(TdA);
    B0xn_gpu = gpuArray(full(B0xn_coeffs));
    B0_gpu = gpuArray(B0);
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
    [A1_gpu,A2_gpu,C1_gpu,C2_gpu,C3_gpu,N_gpu,...
    l1_gpu,l2_gpu,l3_gpu,l4_gpu,l5_gpu,l6_gpu,...
    r1_gpu,r2r3_gpu,...
    shapeDerivative_gpu, l2vec_gpu ]= feval(kernel,...
    RWG.ndof,RWG.ndof,bndmesh.nelt,bndmesh.nvtx,bndmesh.nelt^2,...
    Nthreads,Ivec_gpu,Jvec_gpu,relation_gpu,...
    W0_gpu,X0_gpu,size(X0,1),...
    W1_gpu,X1_gpu,size(X1,1),...
    W2_gpu,X2_gpu,size(X2,1),...
    W3_gpu,X3_gpu,size(X3,1),...
    A1_gpu,A2_gpu,C1_gpu,C2_gpu,C3_gpu,N_gpu,...
    l1_gpu,l2_gpu,l3_gpu,l4_gpu,l5_gpu,l6_gpu,...
    r1_gpu,r2r3_gpu,...
    shapeDerivative_gpu,l2vec_gpu,...
    mu0,mu,...
    TdA_gpu, TnA_gpu,B0xn_gpu,...
    B0_gpu,...
    Elements_gpu,Vertices_gpu,Normals_gpu,Areas_gpu,...
    elt2dof_gpu,elt2dof_gpu,...
    TrialSpace_gpu,TestSpace_gpu,TrialOperator_gpu,TestOperator_gpu,...
    size(RWG.rsf,1),size(RWG.rsf,1),...
    abc_alpha_gpu,Nfields,...
    permI_gpu,permJ_gpu); 


    %% Non SS Computations 

    sdnongpu = zeros(Nfields,1);

    for fieldID = 1:Nfields
        a = abc_alpha(fieldID,1);
        b = abc_alpha(fieldID,2);
        c = abc_alpha(fieldID,3);
        alpha = abc_alpha(fieldID,4);
        [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);

        % Cross bilinear forms
        gradxGdotVelx = @(x,y) 1/4/pi ./vecnorm(x-y,2,2).^3 .* dot(y-x,Vel(x),2);
    
        % ei partial derivative computation
        Aei1mat = integral(Gamma_i,Gamma_e,RWG_i,gradxGdotVelx,RWG_e);
        % int_{Gamma_i} int_{\Gamma_e} gradxG(x,y).vel(x) psi_e(y).zeta_i(x) 
        Aei1 = Psi_i' * Aei1mat * Psi_e;
    
        [Y_e,WY_e] = Gamma_e.qud;
        [X_i,WX_i] = Gamma_i.qud;
        NX = size(X_i,1);
        NY = size(Y_e,1);
        XX = repmat(X_i,NY,1); WWX = repmat(WX_i,NY,1); 
        YY = repelem(Y_e,NX,1); WWY = repelem(WY_e,NX,1);
        W = WWX .* WWY;
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
        % Ciekernel1 = 3/4/pi * (XX-YY).*dot(XX-YY,Vel(XX)-Vel(YY),2)./vecnorm(XX-YY,2,2).^5 ...
        %             -1/4/pi * (Vel(XX)-Vel(YY))./vecnorm(XX-YY,2,2).^3;
    
        % Manually putting vel(x) = 0
        Ciekernel1 = 3/4/pi * (XX-YY).*dot(XX-YY,-Vel(YY),2)./vecnorm(XX-YY,2,2).^5 ...
                    -1/4/pi * (-Vel(YY))./vecnorm(XX-YY,2,2).^3;
    
        Cie1 = sum(W.*dot(Psi_eXX,cross(Ciekernel1,nxgvals_iYY,2),2) ,1);
    
        % 2nd term
        DVel1YY = DVel{1}(YY);
        DVel2YY = DVel{2}(YY);
        DVel3YY = DVel{3}(YY);
    
        DVelnxgvals_iYY = [dot(DVel1YY,nxgvals_iYY,2) dot(DVel2YY,nxgvals_iYY,2) dot(DVel3YY,nxgvals_iYY,2)];
        % Gradx G
        Ciekernel2 = 1/4/pi * (YY-XX)./vecnorm(YY-XX,2,2).^3;
        Cie2 = sum(W.*dot(Psi_eXX,cross(Ciekernel2,DVelnxgvals_iYY,2) ,2) ,1);
    
        [X_i,W_i] = Gamma_i.qud;
        DVelnxgvals_i = [dot(DVel{1}(X_i),nxgvals_i,2)  dot(DVel{2}(X_i),nxgvals_i,2) dot(DVel{3}(X_i),nxgvals_i,2)];
        l7integral = sum(W_i.* DVelnxgvals_i,1);
        l7 = -mu0 * jumpMuInv/2 * dot(B0,l7integral);
        
        Vel_i = Vel(X_i);
        Veldotn_i = dot(Vel_i,normals_i,2);
        r4integral = sum(W_i.*Veldotn_i,1);
        r4 = -jumpMuInv/2 * norm(B0,2)^2 * r4integral;
        
        sdnongpu(fieldID) = -1/(2*mu0) * ( 2 * (Aei1+Aei2)...
                            + 2 * (Cie1+Cie2))...
             +1/mu0 * (l7)...
             + r4;
    end


    sd = shapeDerivative_gpu + sdnongpu;
    
    %pool.delete();
end