% TnA and TdA are the coeffs of the exterior traces of the vector potential
% Psi and g respectively. TnA is expected to be in RWG space.

function sd = SdBemTPVP_dualnorm(bndmesh,TnA,TdA,J,omega_src,mu0,mu,abc_alpha)
    order = 3;
    % BEM Spaces
    P1 = fem(bndmesh,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0 = nxgrad(P1); 
    RWG = fem(bndmesh,'RWG');
    NED = fem(bndmesh,'NED');
    Gamma = dom(bndmesh,order);

    Nfields = size(abc_alpha,1);

    

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
    ptxFilePath = 'SdBemTPVP_dualnorm_GPU.ptx';
    cuFilePath = 'SdBemTPVP_dualnorm_GPU.cu';

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
    A1_gpu = gpuArray.zeros(1,1);
    A2_gpu = gpuArray.zeros(1,1);
    C1_gpu = gpuArray.zeros(1,1);
    C2_gpu = gpuArray.zeros(1,1);
    C3_gpu = gpuArray.zeros(1,1);
    N_gpu = gpuArray.zeros(1,1);
    shapeDerivative_gpu = gpuArray.zeros(Nfields,1);
    TnA_gpu = gpuArray(full(TnA));
    TdA_gpu = gpuArray(TdA);
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
    [A1_gpu,A2_gpu,C1_gpu,C2_gpu,C3_gpu,N_gpu,shapeDerivative_gpu ]= feval(kernel,...
    RWG.ndof,RWG.ndof,bndmesh.nelt,bndmesh.nvtx,bndmesh.nelt^2,...
    Nthreads,Ivec_gpu,Jvec_gpu,relation_gpu,...
    W0_gpu,X0_gpu,size(X0,1),...
    W1_gpu,X1_gpu,size(X1,1),...
    W2_gpu,X2_gpu,size(X2,1),...
    W3_gpu,X3_gpu,size(X3,1),...
    A1_gpu,A2_gpu,C1_gpu,C2_gpu,C3_gpu,N_gpu,...
    shapeDerivative_gpu,...
    mu0,mu,...
    TdA_gpu, TnA_gpu,...
    Elements_gpu,Vertices_gpu,Normals_gpu,Areas_gpu,...
    elt2dof_gpu,elt2dof_gpu,...
    TrialSpace_gpu,TestSpace_gpu,TrialOperator_gpu,TestOperator_gpu,...
    size(RWG.rsf,1),size(RWG.rsf,1),...
    abc_alpha_gpu,Nfields,...
    permI_gpu,permJ_gpu); 
    
    %% non SS computations

    % Getting quadrature points for Gamma and omega_src
    [X,WX] = Gamma.qud;
    [Y,WY] = omega_src.qud;
    % Get tensor product quadrature rule
    NX = size(X,1);
    NY = size(Y,1);
    XX = repmat(X,NY,1); WWX = repmat(WX,NY,1); 
    YY = repelem(Y,NX,1); WWY = repelem(WY,NX,1);
    W = WWX .* WWY;
    JYY = J(YY(:,1),YY(:,2),YY(:,3));

    % Kernel gradx G
    gradxG = 1/(4*pi) * (YY-XX)./vecnorm(XX-YY,2,2).^3;

    gvals = reconstruct(TdA,Gamma,NED);
    normals = Gamma.qudNrm;
    nxgvals = cross(normals,gvals,2);
    nxgvals_XX = repmat(nxgvals,NY,1);

    l1 = zeros(Nfields,1);
    l22 = zeros(Nfields,1);
    l21 = zeros(Nfields,1);

    for fieldID = 1:Nfields
        a = abc_alpha(fieldID,1);
        b = abc_alpha(fieldID,2);
        c = abc_alpha(fieldID,3);
        alpha = abc_alpha(fieldID,4);
        [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);
        % Partial derivative of l1, computed using superconductor shape
        % derivative implementation
        l1(fieldID) = SuperConductorShapeDerivativeT3(bndmesh,TnA,Vel,omega_src,J)...
            + SuperConductorShapeDerivativeT4(bndmesh,TnA,DVel,omega_src,J);
    
        % Partial derivaive of l2
        VelXX = Vel(XX);
        DVel1XX = DVel{1}(XX);
        DVel2XX = DVel{2}(XX);
        DVel3XX = DVel{3}(XX);
        
    
    
        kernell22 = 3/(4*pi) * (XX-YY).*dot(XX-YY,VelXX,2)./vecnorm(XX-YY,2,2).^5 - 1/(4*pi)*VelXX./vecnorm(XX-YY,2,2).^3;
        l22(fieldID) = sum(W.*dot(cross(kernell22,JYY,2),nxgvals_XX,2),1);
    
        
        % DVel {nxg}
        DVelnxgvalsXX = [dot(DVel1XX,nxgvals_XX,2) dot(DVel2XX,nxgvals_XX,2) dot(DVel3XX,nxgvals_XX,2) ];
        l21(fieldID) = sum(W.* dot(cross(gradxG,JYY,2),DVelnxgvalsXX,2),1);
    end

    

    %% SS computations

%     kernelA1 = @(x,y,z) dot(z,Vel(x) - Vel(y), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
% 
%     kernelA2 = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
%     % z := y-x
%     kernelC1 = @(x,y,z) 1/(4*pi) * z./vecnorm(z,2,2).^3;
% 
%     kernelC3 = @(x,y,z) -3/(4*pi) * z .* dot(z,Vel(y)-Vel(x),2)./vecnorm(z,2,2).^5 + 1/(4*pi)*(Vel(y)-Vel(x))./vecnorm(z,2,2).^3;
% 
%     kernelN = kernelA1;
% 
% 
%     Nelt = bndmesh.nelt;
%     [ii,jj] = meshgrid(1:Nelt,1:Nelt);
%     
%     euler = parcluster('local');
%     euler.NumWorkers = 5;
%     saveProfile(euler);
% 
%     pool = euler.parpool(5);
% 
%     spmd
%         if spmdIndex==1
%             % partial derivative of b_A 
%             A1mat = panel_assembly(bndmesh,kernelA1,RWG,RWG,ii(:),jj(:));
%             A1 = TnA' * A1mat * TnA;
% 
%         elseif spmdIndex==2
%             %g
%             DVelRWG = RWG;
%             DVelRWG.opr = 'Dvel[psi]';
%             A2mat = panel_assembly_shape_derivative(bndmesh,kernelA2,DVelRWG,RWG,ii(:),jj(:),Vel,DVel);
%             A2 = 2* TnA' * A2mat * TnA;
% 
%         elseif spmdIndex==3
%             DVelRWG = RWG;
%             DVelRWG.opr = 'Dvel[psi]';
%             % Partial derivative of b_C
%             C1mat = panel_assembly_shape_derivative(bndmesh,kernelC1,DVelRWG,RWG,ii(:),jj(:),Vel,DVel);
%             C1 = TnA' * C1mat * TdA
%             C2 = TdA' * C1mat * TnA % = C1?
% 
%         elseif spmdIndex==4
%             % C3 (Is this way of evaluation okay?), z:= y-x
%             C3mat = panel_assembly(bndmesh,kernelC3,RWG,RWG,ii(:),jj(:));
%             C3 = TnA' * C3mat * TdA
% 
%         elseif spmdIndex==5
%             % Partial derivative of b_N
%             Nmat = panel_assembly_shape_derivative(bndmesh,kernelN,RWG.div,RWG.div,ii(:),jj(:),Vel,DVel);
%             N = -TdA' * Nmat * TdA;
% 
%         end
%     end
    
%     sd = 1/(2*mu0) * ( (1+mu/mu0) * (A1{1}+A2{2})...
%                         -4 * (C1{3}+C2{3}+C3{4})...
%                         +(1+mu0/mu) * (N{5}) )...
%                         -l1 + (l21+l22);

    % sd = 1/(2*mu0) * ( (1+mu/mu0) * (A1+A2)...
    %                     -4 * (C1+C2+C3)...
    %                     +(1+mu0/mu) * (N) )...
    %                     -l1 + (l21+l22);
%     pool.delete();

    sd = 1/(2*mu0) * shapeDerivative_gpu -l1 + (l21+l22);

end
