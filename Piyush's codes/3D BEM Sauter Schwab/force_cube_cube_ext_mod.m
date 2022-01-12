%function [] = force_cube_cube_ext(N)

%close all; clc; clear;
addpath(genpath("../../"));

global X;
global W;

load('X3','X');
load('W3','W');

% Initializing parameters for the problem
N = 40;
%N  = getenv('TESTVAR');
%disp(N);
% Dimensions of the cube
L = [1,1,1];
%distances = [distances; offset];

% Mesh for the geometry
cube_vol_mesh = mshCube(N,L);
mesh_in = cube_vol_mesh.bnd;

cube_vol_mesh = mshCube(N,L);
mesh_out = cube_vol_mesh.bnd;

% Translate the outer cube
tx = 1;
ty = 3;
tz = 2;
dist = norm([tx ty tz]);
N_vtcs = size(mesh_out.vtx,1);
trans = [tx * ones(N_vtcs,1), ty * ones(N_vtcs,1), tz * ones(N_vtcs,1)];
mesh_out.vtx = mesh_out.vtx + trans;

% Join to create the final mesh
mesh = union(mesh_in,mesh_out);

figure;
plot(mesh);
title('Mesh');

%% BEM 

% Definition of FEM spaces and integration rule
S1_Gamma = fem(mesh,'P1');
S0_Gamma = fem(mesh,'P0');
order = 3;
Gamma = dom(mesh,order);
Gamma_in = dom(mesh_in,order);
S0_Gamma_in = fem(mesh_in,'P0');

Op_in = restriction(S0_Gamma,mesh_in);
Op_out = restriction(S0_Gamma,mesh_out);

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
R = 1.1; % Cutoff radius, also used for the velocity field
g = @(X) (sqrt(sum(X.^2,2)) > R)*(2) + (sqrt(sum(X.^2,2)) <= R)*4;

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
Psi_in = Op_in * Psi;
Psi_out = Op_out * Psi;

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
%%
% Classical formula for force evaluation
normals_in = mesh_in.nrm;
classical_force_in = 0.5 * sum( mesh_in.ndv .* Psi_in.^2.* normals_in, 1)

normals_out = mesh_out.nrm;

classical_force_out = 0.5 * sum( mesh_out.ndv .* Psi_out.^2.* normals_out, 1)

charge_in = sum( mesh_in.ndv .* Psi_in, 1)
charge_out = sum( mesh_out.ndv .* Psi_out, 1)
coulomb_force = charge_in * charge_out / 4 / pi / dist^2

%% Computing the classical formula in a different way

% Gives a  1X3 cell
mat = integral(Gamma_in,S0_Gamma_in,ntimes(S0_Gamma_in));

alt_force_in = 0.5*[Psi_in' *mat{1} * Psi_in; Psi_in' *mat{2} * Psi_in; Psi_in' *mat{3} * Psi_in]
