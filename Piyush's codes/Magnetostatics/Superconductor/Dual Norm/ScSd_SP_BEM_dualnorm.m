% Compute the shape derivative obtained from the volume based scalar
% potential formulation of the superconductor

% Super conductor Shape derivative_ Scalar Potential
function sd = ScSd_SP_BEM_dualnorm(bndmesh,Tdu,Tnu,J,omega_src,abc_alpha)
    Nfields = size(abc_alpha,1);
    % BEM Spaces
    P0 = fem(bndmesh,'P0');
    P1 = fem(bndmesh,'P1');
    gradP1 = P1.grad;
    nxgradP1 = P1.nxgrad;

    Gamma = dom(bndmesh,3);

    [X,W] = Gamma.qud;
    normals = Gamma.qudNrm;

    % Reconstructing the Neumann trace psi
    psi = reconstruct(Tnu,Gamma,P0);
    % Reconstructing the surface gradient of u
    gradTu = reconstruct(Tdu,Gamma,gradP1);
    % Reconstructing the trace of u- > g
    g = reconstruct(Tdu,Gamma,P1);

    HJ = compute_vecpot_curl(J,omega_src,X);
    DHJ = compute_vecpot_D_curl(J,omega_src,X);

    %% Evaluating the three terms in the shape derivative

    % Evaluating the velocity field at quadrature points
%     Vels = Vel(X);
%     % Evaluating the Jacobian of the velocity field (row-wise) at qud pts
%     DVel1 = DVel{1}(X);
%     DVel2 = DVel{2}(X);
%     DVel3 = DVel{3}(X);
%     divVel = DVel1(:,1) + DVel2(:,2) + DVel3(:,3);

%     % Stress tensor times normal: T(u) n
%     Tun = 0.5 * (psi.^2 - dot(gradTu,gradTu,2)) .* normals + psi.* gradTu;
%     t1 = sum(W.* dot(Vels, Tun,2) ,1);
% 
%     % 3rd term
%     Vn = dot(Vels,normals,2);
%     t3 = -0.5 * sum(W.* dot(HJ,HJ,2) .* Vn,1);
% 
%     % 2nd term
%     DHJV = [dot(DHJ{1},Vels,2) dot(DHJ{2},Vels,2) dot(DHJ{3},Vels,2)];
%     DVHJ = [dot(DVel1,HJ,2) dot(DVel2,HJ,2) dot(DVel3,HJ,2)];
%     t2 = -sum(W.*g.*dot(normals,(DHJV - DVHJ + HJ.*divVel),2),1);
% 
%     sd = t1+t2+t3;

    %% Full shape derivative
%     Nelt = bndmesh.nelt;

    %% Non panel assembly computations
    
    elseterms = zeros(Nfields,1);

    for fieldID = 1:Nfields
        a = abc_alpha(fieldID,1);
        b = abc_alpha(fieldID,2);
        c = abc_alpha(fieldID,3);
        alpha = abc_alpha(fieldID,4);
        [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);
        
        % t3
        % Evaluating the velocity field at quadrature points
        Vels = Vel(X);
        % Evaluating the Jacobian of the velocity field (row-wise) at qud pts
        DVel1 = DVel{1}(X);
        DVel2 = DVel{2}(X);
        DVel3 = DVel{3}(X);
        divVel = DVel1(:,1) + DVel2(:,2) + DVel3(:,3);
        
        DHJV = [dot(DHJ{1},Vels,2) dot(DHJ{2},Vels,2) dot(DHJ{3},Vels,2)];
        DVHJ = [dot(DVel1,HJ,2) dot(DVel2,HJ,2) dot(DVel3,HJ,2)];
        t3 = -0.5*sum(W.*g.*dot(normals,(DHJV - DVHJ + HJ.*divVel),2),1);
    
        % t4 
        % Projecting the complex integrand to P0 space
        compl_integrand = dot(normals,(DHJV - DVHJ + HJ.*divVel),2);
        compl_integrand_coeffs = proj(compl_integrand,Gamma,P0);
        Kmat = double_layer_laplace(Gamma,P0,P1);
        t4 = -compl_integrand_coeffs' * Kmat * Tdu;
    
        % t8
        HJn = dot(normals,HJ,2);
        HJn_coeffs = proj(HJn,Gamma,P0);
        Vmat = single_layer(Gamma,P0,P0);
        t8 = -HJn_coeffs' * Vmat * compl_integrand_coeffs;
    
        % t9
        Vn = dot(Vels,normals,2);
        t9 = -0.5 * sum(W.* dot(HJ,HJ,2) .* Vn,1);
    
        % t5 and t6
        % Double layer like term in t5
        divVelg = divVel.*g;
        divVelg_coeffs = proj(divVelg,Gamma,P1);
        t5dl = -HJn_coeffs' *Kmat * divVelg_coeffs;
        
        elseterms(fieldID) = t3+t4+t5dl +t8+t9;
    end
    

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
    Ivec = cast(Ivec,'int32');
    Jvec = [J0;J1;J2;J3]-1;
    Jvec = cast(Jvec,'int32');
    relation = [0*ones(size(I0,1),1);...
                1*ones(size(I1,1),1);...
                2*ones(size(I2,1),1);...
                3*ones(size(I3,1),1)];
    relation = cast(relation,'int32');

    % PTX and CUDA files
    ptxFilePath = 'ScSd_SP_BEM_dualnorm_GPU.ptx';
    cuFilePath = 'ScSd_SP_BEM_dualnorm_GPU.cu';

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
    t7_gpu = gpuArray.zeros(1,1);
    t1nt2_gpu = gpuArray.zeros(1,1);
    t56nt6_gpu = gpuArray.zeros(1,1); % excludes the Kmat part which is computed by the CPU
    shapeDerivative_gpu = gpuArray.zeros(Nfields,1);

    % More input variables to GPU
    Tdu_gpu = gpuArray(Tdu);
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
    [t7_gpu,t1nt2_gpu,t56nt6_gpu,shapeDerivative_gpu]= feval(kernel,...
    0,0,bndmesh.nelt,bndmesh.nvtx,bndmesh.nelt^2,...
    Nthreads,Ivec_gpu,Jvec_gpu,relation_gpu,...
    W0_gpu,X0_gpu,size(X0,1),...
    W1_gpu,X1_gpu,size(X1,1),...
    W2_gpu,X2_gpu,size(X2,1),...
    W3_gpu,X3_gpu,size(X3,1),...
    t7_gpu,t1nt2_gpu,t56nt6_gpu,shapeDerivative_gpu,...
    Tdu_gpu, HJn_coeffs_gpu,...
    Elements_gpu,Vertices_gpu,Normals_gpu,Areas_gpu,...
    elt2dof_gpu,elt2dof_gpu,...
    TrialSpace_gpu,TestSpace_gpu,TrialOperator_gpu,TestOperator_gpu,...
    size(P1.rsf,1),size(P1.rsf,1),...
    abc_alpha_gpu,Nfields); 

    %% Panel assembly computations
    % Kernel for t1, z:= y-x
%     [ii,jj] = meshgrid(1:Nelt,1:Nelt);
%     kernelt1 = @(x,y,z) sum(z.*(Vel(x) - Vel(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
%     KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
%     kernelt6 = @(x,y,z) 3/(4*pi)* dot(z,Vel(y) - Vel(x),2) .*z ./vecnorm(z,2,2).^5;
%     combkernel = @(x,y,z) 1/(4*pi) * ( -[ dot(DVel{1}(y),z,2) dot(DVel{2}(y),z,2) dot(DVel{3}(y),z,2) ] + Vel(y) - Vel(x) )./vecnorm(z,2,2).^3;
%     
%     euler = parcluster('local');
%     euler.NumWorkers = 5;
%     saveProfile(euler);
% 
%     pool = euler.parpool(5);

%     spmd
%         if spmdIndex==1
%             % t1
%             t1mat = panel_assembly(bndmesh,kernelt1,nxgradP1,nxgradP1,ii(:),jj(:));
%             t1 = 0.5 * Tdu' * t1mat * Tdu;
% 
%         elseif spmdIndex==2
%             % t2
%             t2mat = panel_assembly_shape_derivative(bndmesh,KV,nxgradP1,nxgradP1,ii(:),jj(:),Vel,DVel);
%             t2 = Tdu' * t2mat * Tdu;
% 
%         elseif spmdIndex==3
%             % t7
%             % Kernel for t1, z:= y-x
%             t7mat = panel_assembly(bndmesh,kernelt1,P0,P0,ii(:),jj(:));
%             t7 = -0.5 * HJn_coeffs' * t7mat * HJn_coeffs;
% 
%         elseif spmdIndex==4
%             % integrable kernel in t6
%             kernelt6mat = panel_assembly(bndmesh,kernelt6,ntimes(P1),P0,ii(:),jj(:));
%             t6 = -HJn_coeffs' * kernelt6mat * Tdu;
% 
%         elseif spmdIndex==5
%             % Combination kernel of t5 and t6 that cancels singularity
%             combkernelmat = panel_assembly(bndmesh,combkernel,ntimes(P1),P0,ii(:),jj(:));
%             t56 = HJn_coeffs' * combkernelmat * Tdu;
%         end
% 
%     end

%     sd_gpu = t1{1}+t2{2}+t56{5}+t6{4}+t7{3};
        sd_gpu_full = shapeDerivative_gpu + elseterms;
%     sd = t1{1}+t2{2}+t3+t4+t5dl+t56{5}+t6{4}+t7{3}+t8+t9;

%     pool.delete();

end