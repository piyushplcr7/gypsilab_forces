function sd_comb = PermanentMagnetShapeDerivativeBIESP_torque(Gamma,Tdu,Tnu,J,omega_src,mu0,M,Xcg)
    % BEM Spaces
    bndmesh = Gamma.msh;
    P0 = fem(bndmesh,'P0');
    P1 = fem(bndmesh,'P1');
    gradP1 = P1.grad;
    nxgradP1 = P1.nxgrad;

    Nfields = 3;

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

    

    %% bilinear form
%     Nelt = bndmesh.nelt;
% 
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

    Mvals = M(Xgypsi);
    Mdotn = dot(Mvals,normals,2);

    lambda_coeffs = proj(Mdotn,Gamma,P0);
    
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
    ptxFilePath = 'PermanentMagnetShapeDerivativeBIESP_torque_GPU.ptx';
    cuFilePath = 'PermanentMagnetShapeDerivativeBIESP_torque_GPU.cu';

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
    l1_gpu = gpuArray.zeros(1,1);
    l2_gpu = gpuArray.zeros(1,1);
    r1_gpu = gpuArray.zeros(1,1);
    shapeDerivative_gpu = gpuArray.zeros(Nfields,1);

    % More input variables to GPU
    Tdu_gpu = gpuArray(Tdu);
    Tnu_gpu = gpuArray(Tnu);
    lambda_coeffs_gpu = gpuArray(full(lambda_coeffs));
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
    Xcg_gpu = gpuArray(Xcg);

    % Launching kernel
    [dbv_ds_gpu,dbw_ds_gpu,dbk_ds_reduced_gpu,l1_gpu,l2_gpu,r1_gpu,shapeDerivative_gpu]= feval(kernel,...
    0,0,bndmesh.nelt,bndmesh.nvtx,bndmesh.nelt^2,...
    Nthreads,Ivec_gpu,Jvec_gpu,relation_gpu,...
    W0_gpu,X0_gpu,size(X0,1),...
    W1_gpu,X1_gpu,size(X1,1),...
    W2_gpu,X2_gpu,size(X2,1),...
    W3_gpu,X3_gpu,size(X3,1),...
    dbv_ds_gpu,dbw_ds_gpu,dbk_ds_reduced_gpu,l1_gpu,l2_gpu,r1_gpu,...
    shapeDerivative_gpu,...
    Tdu_gpu, Tnu_gpu,lambda_coeffs_gpu,mu0,...
    Elements_gpu,Vertices_gpu,Normals_gpu,Areas_gpu,...
    elt2dof_gpu,elt2dof_gpu,...
    TrialSpace_gpu,TestSpace_gpu,TrialOperator_gpu,TestOperator_gpu,...
    size(P1.rsf,1),size(P1.rsf,1),...
    permI_gpu,permJ_gpu,...
    Xcg_gpu); 


    %% CPU

%     euler = parcluster('local');
%     euler.NumWorkers = 5;
%     saveProfile(euler);
% 
%     pool = euler.parpool(5);
% 
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
%             SL_Dvelnxgrad_nxgrad = panel_assembly_shape_derivative(bndmesh,KV,nxgradP1,nxgradP1,ii(:),jj(:),Vel,DVel);
%         end
%     end
    
    Kmat = double_layer_laplace(Gamma,P0,P1);

    % Remaining part of shape derivative computed using cpu
    remaining_parts = zeros(Nfields,1);

    for fieldID = 1:Nfields
        [Vel,DVel] = getRotVelDVel(get_ek(fieldID),Xcg);
        % Evaluating the velocity field at quadrature points
        Vels = Vel(Xgypsi);
        % Evaluating the Jacobian of the velocity field (row-wise) at qud pts
        DVel1 = DVel{1}(Xgypsi);
        DVel2 = DVel{2}(Xgypsi);
        DVel3 = DVel{3}(Xgypsi);
        divVel = DVel1(:,1) + DVel2(:,2) + DVel3(:,3);
        divVelg = divVel.*g;
        divVelg_coeffs = proj(divVelg,Gamma,P1);

        % partial derivative of bk(g,psi), non gpu part
        r2 = mu0 * sum(Wgypsi.* Mdotn .* dot(Vels,HJ,2),1);
        remaining_parts(fieldID) = 2 * mu0 * Tnu' * Kmat * divVelg_coeffs...
            + r2 - mu0 * lambda_coeffs' * Kmat * divVelg_coeffs;
    end
    
    sd_comb = shapeDerivative_gpu + remaining_parts;
    
    
    
end
