% Script for comparing panel_assembly and panel_assembly_vectorized
clear;clc;

%bndmesh = mshSphere(2,1);
mesh = mshCube(2,[1 1 1]);
%mesh = mesh.sub(1);
bndmesh = mesh.bnd;

bndmesh_copy = bndmesh;

P1 = fem(bndmesh,'P1');
P0 = fem(bndmesh,'P0');
% space = space.nxgrad;

elts = 1:bndmesh.nelt;
elts = elts';
elts = repelem(elts,3);
vtcs = reshape(bndmesh.elt',[bndmesh.nelt*3 1]);

Eltmat = sparse(elts,vtcs,ones(size(vtcs)),bndmesh.nelt,bndmesh.nvtx);

Intmat = Eltmat * Eltmat';

Nelt = bndmesh.nelt;

[ii,jj] = meshgrid(1:Nelt,1:Nelt);

KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
KK = @(x,y,z) z./vecnorm(z,2,2).^3/4./pi;

Vnew = panel_assembly_VECTORIZED(bndmesh_copy,KV,P1,P1,Intmat);

Vold = panel_assembly(bndmesh,KV,P1,P1,ii(:),jj(:));

Knew = panel_assembly_VECTORIZED(bndmesh,KK,ntimes(P1),P0,Intmat);

Kold = panel_assembly(bndmesh,KK,ntimes(P1),P0,ii(:),jj(:));



