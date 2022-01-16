close all; clc; clear;
addpath(genpath("../../../../"));
format long;
global X;
global W;
load('X3','X');
load('W3','W');

% Initializing parameters for the problem
N = 100;

% Radius of the sphere
Rad = 5;

% Dimensions of the cube
L = [5 5 5];
for N = 100:50:1700
N
% Mesh for the geometry
mesh_sph = mshSphere(N,Rad);
mesh_cube = mshCube(N,L);
mesh_cube = mesh_cube.bnd;

% Translate the Cube
tx = 10;
ty = 10;
tz = 10;
N_vtcs = size(mesh_cube.vtx,1);
trans = ones(N_vtcs,1) * [tx ty tz];
mesh_cube.vtx = mesh_cube.vtx + trans;

% Translate the sphere
%sph_vtcs = mesh_sph.vtx;
%sph_vtcs(:,1) = sph_vtcs(:,1)*0.5;
%mesh_sph.vtx = sph_vtcs;
mesh_sph.vtx = mesh_sph.vtx + ones(size(mesh_sph.vtx,1),1) * [1 0 0];


% Join to create the final mesh
mesh = union(mesh_cube,mesh_sph);

figure;
plot(mesh);
title('Mesh and normals');
hold on;
% Checking the normal direction for the mesh
ctrs = mesh.ctr;
nrms = mesh.nrm;
quiver3(ctrs(:,1),ctrs(:,2),ctrs(:,3),nrms(:,1),nrms(:,2),nrms(:,3));
scatter3(0,0,0);

%% BEM 

% Definition of FEM spaces and integration rule
S1_Gamma = fem(mesh,'P1');
S0_Gamma = fem(mesh,'P0');
order = 3;
Gamma = dom(mesh,order);

Proj_Op_cube = restriction(S0_Gamma,mesh_cube);
Proj_Op_sph = restriction(S0_Gamma,mesh_sph);

% Solving a Direct first kind BVP to get the representation of the state
% V Psi = (0.5*M+K) g_N (Interior problem)
% V Psi = (-0.5*M+K) g_N (Exterior problem)

% Getting the Single Layer matrix V using Gypsilab implementation
Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); % 0 wave number
V = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,Gxy,S0_Gamma);
V = V + 1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'[1/r]',S0_Gamma);

% Getting the Double Layer matrix K using Gypsilab implementation
GradG = cell(3,1);
GradG{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
GradG{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
GradG{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);
K = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,GradG,ntimes(S0_Gamma));
K = K +1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'grady[1/r]',ntimes(S0_Gamma));

% Defining the mass matrix M
M = integral(Gamma,S0_Gamma,S0_Gamma);

% Defining the Dirichlet boundary condition
g = @(x) 1/4/pi./vecnorm(x,2,2); % Point charge at origin

figure;
plot(mesh);
hold on;
plot(mesh,g(S0_Gamma.dof));
title('Dirichlet Boundary Condition');
colorbar;

% Constructing the RHS
g_N = integral(Gamma,S0_Gamma,g);

% Exterior problem
Psi = V\((-0.5 * M + K)* (M\g_N));

% Getting Neumann trace on the two distinct surfaces
Psi_cube = Proj_Op_cube * Psi;
Psi_sph = Proj_Op_sph * Psi;

% Solving the adjoint problem to get the adjoint solution
Rho = V\(-0.5 * g_N);

% Visualizing the Neumann trace
dofs = S0_Gamma.dof;
normals = mesh.nrm;
sPsi = -2*Psi;
figure;
plot(mesh);
hold on;
quiver3(dofs(:,1),dofs(:,2),dofs(:,3),normals(:,1).*sPsi,normals(:,2).*sPsi,normals(:,3).*sPsi,0);
title('Field on the surface');
%quiver3(dofs(:,1),dofs(:,2),dofs(:,3),normals(:,1),normals(:,2),normals(:,3));

% Visualizing the charge density on the surface of the objects
figure;
plot(mesh);
hold on;
plot(mesh,Psi);
title("Surface charge density");
colorbar;

%% Exact solution
Psi_exact = -1/4/pi./vecnorm(dofs,2,2).^3 .* (dot(dofs,nrms,2));
Psi_cube_exact = Proj_Op_cube * Psi_exact;
Psi_sph_exact = Proj_Op_sph * Psi_exact;

% Visualizing the exactness
figure;
plot(Psi./Psi_exact);
title('Psi/Psi exact');

figure;
plot(Psi_sph./Psi_sph_exact);
title('Psi_sph/Psi sph exact');

figure;
plot(Psi_cube./Psi_cube_exact);
title('Psi tor/Psi tor exact');

err_vec = Psi-Psi_exact;

% L2 Error
l2_err = err_vec' * M * err_vec%norm(err_vec)

av_err = err_vec' * V * err_vec

linf_err = max(abs(err_vec))

close all;
end