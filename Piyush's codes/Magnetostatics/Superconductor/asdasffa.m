clear; clc;

i = 1;
vals = [1 1/2 1/4 1/7.9 1/15.9];
N = 2^(i+6); 
N_src = N;
R0 = 2;
r0 = .5;
[J,mesh_src] = get_torus_source(N_src,R0,r0);
omega_src = dom(mesh_src,3);

% Old tetra mesh
bndmesh = meshSymTetra;
bndmesh = bndmesh.translate([2 1 3]);
bndmesh = bndmesh.refine(vals(i));
col = "red";


% % New tetra mesh
tetra_function_name = sprintf('tetra%d', i);
tetra_function_handle = str2func(tetra_function_name);
bndmesh = genMeshFromScript(tetra_function_handle);
bndmesh = bndmesh.translate([2 1 3]);
col = "blue";

% Rest of the procedure
vtcs = bndmesh.vtx;
Gamma = dom(bndmesh,3);
[X,~,~] = Gamma.qud;
[Tdu,Tnu] = solve_superconductor_scalar_potential(bndmesh,J,mesh_src);

P1 = fem(bndmesh,'P1');
SCurlP1 = nxgrad(P1);

normals = Gamma.qudNrm;

quiver3wrapper(X,reconstruct(Tdu,Gamma,SCurlP1),col); hold on;

% scatter3(X(:,1),X(:,2),X(:,3),col); hold on;

% scatter3(vtcs(:,1),vtcs(:,2),vtcs(:,3),col); hold on;

% quiver3wrapper(X,normals,col); hold on;

% testmesh = mshTetra(1);