% Script testing CUDA
addpath(genpath("../../../"));

clear;clc;

format long;

%bndmesh = mshSphere(2,1);
mesh = mshCube(2000,[1 1 1]);
% mesh = mesh.sub(1);
bndmesh = mesh.bnd;

% bndmesh_copy = bndmesh;

P1 = fem(bndmesh,'P1');
P0 = fem(bndmesh,'P0');

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

ptxFilePath = '../../CUDA/SauterSchwabQuadrature.ptx';
cuFilePath = '../../CUDA/SauterSchwabQuadrature.cu';
kernelName = 'computeShapeDerivative';

order = 10;
[X, W] = quad4D(order);

% CUDAKernel object
% kernel = parallel.gpu.CUDAKernel(ptxFilePath, kernelName);
% kernel = parallel.gpu.CUDAKernel(ptxFilePath, cuFilePath,kernelName);
kernel = parallel.gpu.CUDAKernel(ptxFilePath, cuFilePath);

% Set grid and block dimensions (example values).
gridDim = [100, 1, 1];
blockDim = [512, 1, 1];

% gridDim = [1, 1, 1];
% blockDim = [1, 1, 1];

Nthreads = prod(blockDim) * prod(gridDim);

kernel.GridSize = gridDim;
kernel.ThreadBlockSize = blockDim;

% Transferring inputs to the GPU
Ivec_gpu = gpuArray(Ivec);
Jvec_gpu = gpuArray(Jvec);
relation_gpu = gpuArray(relation);

W0_gpu = gpuArray(W{4});
W1_gpu = gpuArray(W{3});
W2_gpu = gpuArray(W{2});
W3_gpu = gpuArray(W{1});

X0 = X{4}; X0(:,1) = X0(:,1)-X0(:,2); X0(:,3) = X0(:,3)-X0(:,4);
X0_gpu = gpuArray(X0);

X1 = X{3}; X1(:,1) = X1(:,1)-X1(:,2); X1(:,3) = X1(:,3)-X1(:,4);
X1_gpu = gpuArray(X1);

X2 = X{2}; X2(:,1) = X2(:,1)-X2(:,2); X2(:,3) = X2(:,3)-X2(:,4);
X2_gpu = gpuArray(X2);

X3 = X{1}; X3(:,1) = X3(:,1)-X3(:,2); X3(:,3) = X3(:,3)-X3(:,4);
X3_gpu = gpuArray(X3);

shapeDerivative_gpu = gpuArray.zeros(1,1);

GalerkinMatrix_gpu = gpuArray.zeros(P0.ndof,P0.ndof);

trialVec_gpu = gpuArray.zeros(1,1);

testVec_gpu = gpuArray.zeros(1,1);

Elements = cast(bndmesh.elt-1,'int32');

Elements_gpu = gpuArray(Elements);

Vertices_gpu = gpuArray(bndmesh.vtx);

Normals_gpu = gpuArray(bndmesh.nrm);

Areas_gpu = gpuArray(bndmesh.ndv);

TrialSpace_gpu = 0;

TestSpace_gpu = 0;

TrialOperator_gpu = 0;

TestOperator_gpu = 0;

NRSFTrial = size(P0.rsf,1);

NRSFTest = size(P0.rsf,1);

%%

% testout_gpu = gpuArray.zeros(size(W{2},1),1);
% testABCi_gpu = gpuArray.zeros(3,3); 
% testABCj_gpu = gpuArray.zeros(3,3); 
% origelti_gpu = gpuArray.zeros(3,1,'int32');
% origeltj_gpu = gpuArray.zeros(3,1,'int32');
% modifelti_gpu = gpuArray.zeros(3,1,'int32');
% modifeltj_gpu = gpuArray.zeros(3,1,'int32');

% Launch CUDA kernel
[sd,gmat] = feval(kernel,...
    P0.ndof,P0.ndof,bndmesh.nelt,bndmesh.nvtx,...
    Nthreads,Ivec_gpu,Jvec_gpu,relation_gpu,...
    W0_gpu,X0_gpu,size(X0,1),...
    W1_gpu,X1_gpu,size(X1,1),...
    W2_gpu,X2_gpu,size(X2,1),...
    W3_gpu,X3_gpu,size(X3,1),...
    shapeDerivative_gpu,...
    GalerkinMatrix_gpu,...
    trialVec_gpu, testVec_gpu,...
    Elements_gpu,Vertices_gpu,Normals_gpu,Areas_gpu,...
    TrialSpace_gpu,TestSpace_gpu,TrialOperator_gpu,TestOperator_gpu,...
    NRSFTrial,NRSFTest); %...
%     ,testout_gpu,testABCi_gpu,testABCj_gpu, ...
%     origelti_gpu,origeltj_gpu,modifelti_gpu,modifeltj_gpu);

KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
Nelt = bndmesh.nelt;
[ii,jj] = meshgrid(1:Nelt,1:Nelt);

tic;
V = panel_assembly(bndmesh,KV,P0,P0,ii(:),jj(:));

elapsed_time = toc;
fprintf('CPU took %.4f seconds.\n', elapsed_time);

