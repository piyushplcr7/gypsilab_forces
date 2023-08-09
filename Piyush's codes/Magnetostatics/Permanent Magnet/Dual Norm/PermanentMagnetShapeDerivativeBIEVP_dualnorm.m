% This function implements the shape derivative assuming a constant
% magnetization

% Inputs
%
% TdA       : Coefficients for the exterior Dirichlet trace, in NED space
% TnA       : Coefficients for the exterior Neumann trace, in RWG

function sd = PermanentMagnetShapeDerivativeBIEVP_dualnorm(Gamma,TnA,TdA,J,omega_src,mu,mu0,M,abc_alpha)
    %% BEM Spaces
    % Source mesh
    % Magnet boundary mesh
    bndmesh = Gamma.msh;
    [X,~] = Gamma.qud;
    normals = Gamma.qudNrm;
     % BEM spaces
    % curl conforming -> Dirichlet trace
    NED = fem(bndmesh,'NED'); 
    % Div conforming space 
    RWG = fem(bndmesh,'RWG');

    Nfields = size(abc_alpha,1);

    %% Remaining terms + partial derivative of l_M
    Nelt = bndmesh.nelt;
    [ii,jj] = meshgrid(1:Nelt,1:Nelt);

    Mvals = M(X);
    Mxn = cross(Mvals,normals,2);
    Mxn_coeffs = proj(Mxn,Gamma,RWG);
    
%     kernelA1 = @(x,y,z) dot(z,Vel(x) - Vel(y), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
%     kernelA2 = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
%     kernelC1 = @(x,y,z) 1/(4*pi) * z./vecnorm(z,2,2).^3;
%     kernelC3 = @(x,y,z) -3/(4*pi) * z .* dot(z,Vel(y)-Vel(x),2)./vecnorm(z,2,2).^5 + 1/(4*pi)*(Vel(y)-Vel(x))./vecnorm(z,2,2).^3;
%     kernelN = kernelA1;
    
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
    ptxFilePath = 'PermanentMagnetShapeDerivativeBIEVP_dualnorm_GPU.ptx';
    cuFilePath = 'PermanentMagnetShapeDerivativeBIEVP_dualnorm_GPU.cu';

    % Quadrature to be passed to the GPU
    order = 5;
    [X, W] = quad4D(order);
    
    % Kernel Object
    kernel = parallel.gpu.CUDAKernel(ptxFilePath, cuFilePath);
    gridDim = [150, 1, 1];
    blockDim = [32, 1, 1]; % 32 is the SIMD Width or wrap size
    Nthreads = prod(blockDim) * prod(gridDim);
    kernel.GridSize = gridDim;
    kernel.ThreadBlockSize = blockDim; 
    
    % Transferring Interactions and Relations to GPU
    Ivec_gpu = gpuArray(Ivec);
    Jvec_gpu = gpuArray(Jvec);
    relation_gpu = gpuArray(relation);
    
    % Transferring Quadrature to GPU
    W0 = W{4};
    W1 = W{3};
    W2 = W{2};
    W3 = W{1};
    W0_gpu = gpuArray(W0);
    W1_gpu = gpuArray(W1);
    W2_gpu = gpuArray(W2);
    W3_gpu = gpuArray(W3);
    X0 = X{4}; X0(:,1) = X0(:,1)-X0(:,2); X0(:,3) = X0(:,3)-X0(:,4);
    X0 = X0';
    X0_gpu = gpuArray(X0);
    X1 = X{3}; X1(:,1) = X1(:,1)-X1(:,2); X1(:,3) = X1(:,3)-X1(:,4);
    X1 = X1';
    X1_gpu = gpuArray(X1);
    X2 = X{2}; X2(:,1) = X2(:,1)-X2(:,2); X2(:,3) = X2(:,3)-X2(:,4);
    X2 = X2';
    X2_gpu = gpuArray(X2);
    X3 = X{1}; X3(:,1) = X3(:,1)-X3(:,2); X3(:,3) = X3(:,3)-X3(:,4);
    X3 = X3';
    X3_gpu = gpuArray(X3);

    % More variables for GPU
    A1_gpu = gpuArray.zeros(1,1);
    A2_gpu = gpuArray.zeros(1,1);
    C1_gpu = gpuArray.zeros(1,1);
    C2_gpu = gpuArray.zeros(1,1);
    C3_gpu = gpuArray.zeros(1,1);
    MC1_gpu = gpuArray.zeros(1,1);
    MC2_gpu = gpuArray.zeros(1,1);
    MC3_gpu = gpuArray.zeros(1,1);
    N_gpu = gpuArray.zeros(1,1);
    shapeDerivative_gpu = gpuArray.zeros(Nfields,1);
    red_remaining_gpu = gpuArray.zeros(1,1);
    red_l_M_gpu = gpuArray.zeros(1,1);
    blue_remaining_gpu = gpuArray.zeros(1,1);
    blue_l_M_gpu = gpuArray.zeros(1,1);
    TnA_gpu = gpuArray(full(TnA));
    TdA_gpu = gpuArray(TdA);
    Mxn_coeffs_gpu = gpuArray(full(Mxn_coeffs));
    Elements = cast(bndmesh.elt'-1,'int32');
    Elements_gpu = gpuArray(Elements);
    Vertices = bndmesh.vtx';
    Vertices_gpu = gpuArray(Vertices);
    Normals = bndmesh.nrm';
    Normals_gpu = gpuArray(Normals);
    Areas = bndmesh.ndv;
    Areas_gpu = gpuArray(Areas);
    zeroIdxelt2dof = cast(elt2dof'-1,'int32');
    elt2dof_gpu = gpuArray(zeroIdxelt2dof);
    TrialSpace_gpu = 0;
    TestSpace_gpu = 0;
    TrialOperator_gpu = 0;
    TestOperator_gpu = 0;
    abc_alpha_gpu = gpuArray(cast(abc_alpha','int32'));

    % Launching kernel
    [A1_gpu,A2_gpu,C1_gpu,C2_gpu,C3_gpu,N_gpu,red_remaining_gpu, red_l_M_gpu, blue_remaining_gpu, blue_l_M_gpu,MC1_gpu,MC2_gpu,MC3_gpu,shapeDerivative_gpu ]= feval(kernel,...
    RWG.ndof,RWG.ndof,bndmesh.nelt,bndmesh.nvtx,bndmesh.nelt^2,...
    Nthreads,Ivec_gpu,Jvec_gpu,relation_gpu,...
    W0_gpu,X0_gpu,size(X0,2),...
    W1_gpu,X1_gpu,size(X1,2),...
    W2_gpu,X2_gpu,size(X2,2),...
    W3_gpu,X3_gpu,size(X3,2),...
    A1_gpu,A2_gpu,C1_gpu,C2_gpu,C3_gpu,N_gpu,...
    red_remaining_gpu, red_l_M_gpu, blue_remaining_gpu, blue_l_M_gpu,...
    MC1_gpu,MC2_gpu,MC3_gpu,...
    shapeDerivative_gpu,...
    TdA_gpu, TnA_gpu,Mxn_coeffs_gpu,...
    mu,mu0,...
    Elements_gpu,Vertices_gpu,Normals_gpu,Areas_gpu,...
    elt2dof_gpu,elt2dof_gpu,...
    TrialSpace_gpu,TestSpace_gpu,TrialOperator_gpu,TestOperator_gpu,...
    size(RWG.rsf,1),size(RWG.rsf,1),...
    abc_alpha_gpu,Nfields,...
    permI_gpu,permJ_gpu); 

    %% MEX CPU

    shapeDerivative_cpu = PermanentMagnetShapeDerivativeBIEVP_dualnorm_CPU(...
        cast(RWG.ndof,'int32'),...
        cast(RWG.ndof,'int32'),...
        cast(bndmesh.nelt,'int32'),...
        cast(bndmesh.nvtx,'int32'),...
        cast(size(Ivec,1),'int32'),...
        Ivec,...
        Jvec,...
        relation,...
        W0,...
        X0,...
        cast(size(X0,2),'int32'),...
        W1,...
        X1,...
        cast(size(X1,2),'int32'),...
        W2,...
        X2,...
        cast(size(X2,2),'int32'),...
        W3,...
        X3,...
        cast(size(X3,2),'int32'),...
        0,...
        0,...
        Elements,...
        Vertices,...
        Normals,...
        Areas,...
        zeroIdxelt2dof,...
        zeroIdxelt2dof,...
        cast(1,'int32'),... % Trial Space
        cast(1,'int32'),... % Test Space
        cast(0,'int32'),... % Trial Op
        cast(0,'int32'),... % Test Op
        cast(size(RWG.rsf,1),'int32'),...
        cast(size(RWG.rsf,1),'int32'),...
        TdA,...
        full(TnA),...
        full(Mxn_coeffs),...
        cast(abc_alpha','int32'),...
        cast(Nfields,'int32'),...
        mu,...
        mu0...
        );

    %% CPU
%     DVelRWG = RWG;
%     DVelRWG.opr = 'Dvel[psi]';
% 
%     euler = parcluster('local');
%     euler.NumWorkers = 5;
%     saveProfile(euler);
% 
%     pool = euler.parpool(5);
% 
%     spmd
%         if spmdIndex==1
%             % Red term in implementation.xopp, pg 67
%             % z := y-x, grad_x G = z/norm(z)^3
%             A1mat = panel_assembly(bndmesh,kernelA1,RWG,RWG,ii(:),jj(:));
%             red_remaining = mu/2 * Mxn_coeffs' * A1mat * Mxn_coeffs;
%             red_l_M = mu * TnA' * A1mat * Mxn_coeffs;
%             A1 = TnA' * A1mat * TnA;
% 
%         elseif spmdIndex==2
%             % Blue terms in implementation.xopp, pg 67
%             A2mat = panel_assembly_shape_derivative(bndmesh,kernelA2,DVelRWG,RWG,ii(:),jj(:),Vel,DVel);
%             blue_remaining = mu * Mxn_coeffs' * A2mat * Mxn_coeffs;
%             blue_l_M = mu * (Mxn_coeffs' * A2mat * TnA + TnA' * A2mat * Mxn_coeffs);
%             A2 = 2* TnA' * A2mat * TnA;
% 
%         elseif spmdIndex==3
%             % Derivative of double layer in l_M, borrowed from TP
%             % Partial derivative of b_C
%             % z := y-x
%             C1mat = panel_assembly_shape_derivative(bndmesh,kernelC1,DVelRWG,RWG,ii(:),jj(:),Vel,DVel);
%             MC1 = -mu0* Mxn_coeffs' * C1mat * TdA;
%             MC2 = -mu0* TdA' * C1mat * Mxn_coeffs; 
%             C1 = TnA' * C1mat * TdA;
%             C2 = TdA' * C1mat * TnA; 
% 
%         elseif spmdIndex==4
%             % C3 (Is this way of evaluation okay?), z:= y-x
%             C3mat = panel_assembly(bndmesh,kernelC3,RWG,RWG,ii(:),jj(:));
%             MC3 = -mu0* Mxn_coeffs' * C3mat * TdA;
%             C3 = TnA' * C3mat * TdA;
% 
%         elseif spmdIndex==5
%             % Partial derivative of b_N
%             Nmat = panel_assembly_shape_derivative(bndmesh,kernelN,RWG.div,RWG.div,ii(:),jj(:),Vel,DVel);
%             N = -TdA' * Nmat * TdA;
% 
%         end
% 
%     end

    %% Shape derivative of linear forms from transmission problem

    l1 = zeros(Nfields,1);
    l22 = zeros(Nfields,1);
    l21 = zeros(Nfields,1);

    % Partial derivaive of l2
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
    gvals = reconstruct(TdA,Gamma,NED);
    nxgvals = cross(normals,gvals,2);
    nxgvals_XX = repmat(nxgvals,NY,1);
    % Kernel gradx G
    gradxG = 1/(4*pi) * (YY-XX)./vecnorm(XX-YY,2,2).^3;

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
    
    % Combining the shape derivative from different sources
    sd = shapeDerivative_gpu -l1 + l21 + l22;

    

%     % SD from transmission problem
%     sd = 1/(2*mu0) * ( (1+mu/mu0) * (A1{1}+A2{2})...
%                         -4 * (C1{3}+C2{3}+C3{4})...
%                         +(1+mu0/mu) * (N{5}) )...
%                         -l1 + (l21+l22); 
%     % SD additional terms for magnet
%     sd = sd + red_remaining{1} + red_l_M{1} + blue_remaining{2} + blue_l_M{2} + MC1{3} + MC2{3} + MC3{4};
%     pool.delete();
end