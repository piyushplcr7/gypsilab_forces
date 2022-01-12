%function l2_err = compute_err_ext_sph(N)
close all; clc; clear;
addpath(genpath("../../"));

global X;
global W;
load('X3','X');
load('W3','W');

% Initializing parameters for the problem
% Radius of spheres
Rad = 10;
% distance between the centers
% separation
s = Rad;
dist = 2*Rad + s;

% Mesh for the geometry
mesh_in = mshSphere(N,Rad);
mesh_out = mesh_in;
% Translate the outer cube
N_vtcs = size(mesh_out.vtx,1);
trans = [dist * ones(N_vtcs,1), zeros(N_vtcs,1), zeros(N_vtcs,1)];
mesh_out.vtx = mesh_out.vtx + trans;

% Slightly translating the inner sphere
mesh_in.vtx = mesh_in.vtx + 0.1*trans;

% Join to create the final mesh
mesh = union(mesh_in,mesh_out);

figure;
plot(mesh); 
title('Mesh and normals');
hold on;
% Checking the normal direction for the mesh
ctrs = mesh.ctr;
nrms = mesh.nrm;
quiver3(ctrs(:,1),ctrs(:,2),ctrs(:,3),nrms(:,1),nrms(:,2),nrms(:,3));

%% Solving the Laplace problem using BEM

% Definition of FEM spaces and integration rule
S1_Gamma = fem(mesh,'P1');
S0_Gamma = fem(mesh,'P0');
order = 3;
Gamma = dom(mesh,order);

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
g = @(x) 1/4/pi./vecnorm(x,2,2); % Point charge at origin

% Checking the boundary condition
figure;
plot(mesh);
hold on;
PI_g = g(S0_Gamma.dof);
plot(mesh,PI_g);
title('Boundary condition');
colorbar;

% Constructing the RHS
g_N = integral(Gamma,S0_Gamma,g);

% Exterior problem
Psi = V\((-0.5 * M + K)* (M\g_N));

Psi_in = Op_in * Psi;
Psi_out = Op_out * Psi;

% Visualizing the Neumann trace
dofs = S0_Gamma.dof;
normals = mesh.nrm;
sPsi = 2*Psi;
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
Psi_exact = -1/4/pi./vecnorm(dofs,2,2).^3 .* (dofs.*nrms);

% Visualizing the exactness
figure;
plot(Psi./Psi_exact);

% L2 Error
l2_err = norm(Psi-Psi_exact);
close all;
%end