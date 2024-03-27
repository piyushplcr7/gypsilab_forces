% Script to compare the operator evaluated using gypsilab and mex file
clear;clc;
addpath(genpath("../../../"));
N = 200;
L = 2*[1 1 1];
T = [5 5 3];

mesh = mshCube(N,L);
mesh = mesh.translate(T);
bndmesh = mesh.bnd;

P1 = fem(bndmesh,'P1');
P0 = fem(bndmesh,'P0');

Gamma = dom(bndmesh,3);

V_gypsi = single_layer(Gamma,P0,P0);
K_gypsi = double_layer_laplace(Gamma,P0,P1);
W_gypsi = single_layer(Gamma,P1.nxgrad,P1.nxgrad);

%% Calling the mex function
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
Ivec = cast(Ivec,'int32');
Jvec = [J0;J1;J2;J3]-1;
%     Jvec = cast(Jvec,'int32');
Jvec = cast(Jvec,'int32');
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

permI_gpu = cast(permI','int32');
permJ_gpu = cast(permJ','int32');

order = 5;
[X, W] = quad4D(order);

W0_gpu = W{4};
W1_gpu = W{3};
W2_gpu = W{2};
W3_gpu = W{1};

X0 = X{4}; X0(:,1) = X0(:,1)-X0(:,2); X0(:,3) = X0(:,3)-X0(:,4);
X0_gpu = X0';
X1 = X{3}; X1(:,1) = X1(:,1)-X1(:,2); X1(:,3) = X1(:,3)-X1(:,4);
X1_gpu = X1';
X2 = X{2}; X2(:,1) = X2(:,1)-X2(:,2); X2(:,3) = X2(:,3)-X2(:,4);
X2_gpu = X2';
X3 = X{1}; X3(:,1) = X3(:,1)-X3(:,2); X3(:,3) = X3(:,3)-X3(:,4);
X3_gpu = X3';

Elements = cast(bndmesh.elt-1,'int32');
Elements_gpu = Elements';
Vertices_gpu = bndmesh.vtx';
Normals_gpu = bndmesh.nrm';
Areas_gpu = bndmesh.ndv;

tic;
[V,K,W] = LaplaceTPOperator(cast(bndmesh.nelt,'int32'), cast(bndmesh.nvtx,'int32'), cast(bndmesh.nelt^2,'int32'),... 
    Ivec, Jvec, relation, ...
    W0_gpu,X0_gpu,cast(size(X0,1),'int32'),...
    W1_gpu,X1_gpu,cast(size(X1,1),'int32'),...
    W2_gpu,X2_gpu, cast(size(X2,1),'int32'),...
    W3_gpu,X3_gpu, cast(size(X3,1),'int32'),...
    Elements_gpu,Vertices_gpu,Normals_gpu,Areas_gpu,...
    permI_gpu,permJ_gpu);

time_elapsed_1 = toc

tic;
[VP,KP,WP] = LaplaceTPOperatorParallel(cast(bndmesh.nelt,'int32'), cast(bndmesh.nvtx,'int32'), cast(bndmesh.nelt^2,'int32'),... 
    Ivec, Jvec, relation, ...
    W0_gpu,X0_gpu,cast(size(X0,1),'int32'),...
    W1_gpu,X1_gpu,cast(size(X1,1),'int32'),...
    W2_gpu,X2_gpu, cast(size(X2,1),'int32'),...
    W3_gpu,X3_gpu, cast(size(X3,1),'int32'),...
    Elements_gpu,Vertices_gpu,Normals_gpu,Areas_gpu,...
    permI_gpu,permJ_gpu);

time_elapsed_par = toc



%% Assembling using SS in gypsi
 Nelt = bndmesh.nelt;
[ii,jj] = meshgrid(1:Nelt,1:Nelt);

SLKernel = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
DLKernel = @(x,y,z) -1/4/pi * z./vecnorm(z,2,2).^3;

tic;

VSS_gypsi = panel_assembly(bndmesh,SLKernel,P0,P0,ii(:),jj(:));
KSS_gypsi = panel_assembly(bndmesh,DLKernel,ntimes(P1),P0,ii(:),jj(:));
WSS_gypsi = panel_assembly(bndmesh,SLKernel,P1.nxgrad,P1.nxgrad,ii(:),jj(:));

time_elapsed_2 = toc





