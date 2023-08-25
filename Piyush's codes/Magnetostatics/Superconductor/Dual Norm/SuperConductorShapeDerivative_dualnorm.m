% Superconductor shape derivative for Neumann Trace input lying in RWG
% Space

function val_gpu = SuperConductorShapeDerivative_dualnorm(bndmesh,TnA,omega_src,J,abc_alpha)
    
    RWG = fem(bndmesh,'RWG');
    % Number of fields
    Nfields = size(abc_alpha,1);

    % Kernel, z:= y-x
%     kernelgypsi = @(x,y,z) sum(z.*(Vel(x) - Vel(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);

%     Nelt = bndmesh.nelt;

%     [ii,jj] = meshgrid(1:Nelt,1:Nelt);
    
    %% Dispatching SS computations to GPU

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
    ptxFilePath = 'SuperConductorShapeDerivative_dualnorm_GPU.ptx';
    cuFilePath = 'SuperConductorShapeDerivative_dualnorm_GPU.cu';

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

    % More variables for GPU
    val1_gpu = gpuArray.zeros(1,1);
    val2_gpu = gpuArray.zeros(1,1);
    shapeDerivative_gpu = gpuArray.zeros(Nfields,1);

    TnA_gpu = gpuArray(full(TnA));
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
    [val1_gpu,val2_gpu,shapeDerivative_gpu]= feval(kernel,...
    RWG.ndof,RWG.ndof,bndmesh.nelt,bndmesh.nvtx,bndmesh.nelt^2,...
    Nthreads,Ivec_gpu,Jvec_gpu,relation_gpu,...
    W0_gpu,X0_gpu,size(X0,1),...
    W1_gpu,X1_gpu,size(X1,1),...
    W2_gpu,X2_gpu,size(X2,1),...
    W3_gpu,X3_gpu,size(X3,1),...
    val1_gpu,val2_gpu,shapeDerivative_gpu,...
    TnA_gpu,...
    Elements_gpu,Vertices_gpu,Normals_gpu,Areas_gpu,...
    elt2dof_gpu,elt2dof_gpu,...
    TrialSpace_gpu,TestSpace_gpu,TrialOperator_gpu,TestOperator_gpu,...
    size(RWG.rsf,1),size(RWG.rsf,1),...
    abc_alpha_gpu,Nfields,...
    permI_gpu,permJ_gpu); 

    %% CPU

%     t1mat = panel_assembly(bndmesh,kernelgypsi,RWG,RWG,ii(:),jj(:));
%     % 1st term
%     val = 0.5 * dot(TnA,t1mat*TnA);
% 
%     % 2nd term 
%     % Defining the kernel for single layer BIO
%     % KV = @(x,y,z) 1/norm(z)/4./pi;
%     KV = @(x,y,z) sqrt(1./ sum(z.^2 ,2) ) /4./pi;
% 
%     DVelRWG = RWG;
%     DVelRWG.opr = 'Dvel[psi]';
% 
%     t2mat = panel_assembly_shape_derivative(bndmesh,KV,DVelRWG,RWG,ii(:),jj(:),Vel,DVel);
%     val = val + dot(TnA,t2mat*TnA);

    t3 = zeros(Nfields,1);
    t4 = zeros(Nfields,1);

    for fieldID = 1:Nfields
        a = abc_alpha(fieldID,1);
        b = abc_alpha(fieldID,2);
        c = abc_alpha(fieldID,3);
        alpha = abc_alpha(fieldID,4);
        [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);
        t3(fieldID) = SuperConductorShapeDerivativeT3(bndmesh,TnA,Vel,omega_src,J);
        t4(fieldID) = SuperConductorShapeDerivativeT4(bndmesh,TnA,DVel,omega_src,J);
        
    end

    % 3rd term
%     t3 = SuperConductorShapeDerivativeT3(bndmesh,TnA,Vel,omega_src,J);
% 
%     % 4th term
%     t4 = SuperConductorShapeDerivativeT4(bndmesh,TnA,DVel,omega_src,J);

%     val = val - t3-t4;

    val_gpu = shapeDerivative_gpu -t3 -t4;
end
