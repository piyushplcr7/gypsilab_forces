%function [] = force_cube_cube_ext(N)

close all; clc; clear;
addpath(genpath("../../../"));
format long;
global X;
global W;

load('X3','X');
load('W3','W');

% Initializing parameters for the problem
N = 100;

% Radius of the sphere
Rad = 5;

% Radii for the torus
r1 = 5;
r2 = 2;
for N = 100:100:3000
N
% Mesh for the geometry
mesh_sph = mshSphere(N,Rad);

mesh_tor = mshTorus(N,r1,r2);

% Rotate the torus
rng(10);
A = rand(3,3);
[Q,R] = qr(A);
mesh_tor.vtx = mesh_tor.vtx * Q;

% Translate the torus
tx = 15;
ty = 0;
tz = 0;
N_vtcs = size(mesh_tor.vtx,1);
trans = ones(N_vtcs,1) * [tx ty tz];
mesh_tor.vtx = mesh_tor.vtx + trans;

% Join to create the final mesh
mesh = union(mesh_tor,mesh_sph);

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

Proj_Op_tor = restriction(S0_Gamma,mesh_tor);
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
R = 8.5; % Cutoff radius, also used for the velocity field
g = @(X) (sqrt(sum(X.^2,2)) > R)*(2) + (sqrt(sum(X.^2,2)) <= R)*4;
%g = @(x) 1/4/pi./vecnorm(x,2,2); % Point charge at origin

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
Psi_tor = Proj_Op_tor * Psi;
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
%% Checking the convergence of the Neumann trace at a point
testpt = [Rad 0 0];
dofs = S0_Gamma.dof;
Ndofs = size(dofs,1);
distances = vecnorm((dofs - ones(Ndofs,1) * testpt),2,2);
[min_d,min_ind] = min(distances);
Psi_testpt = Psi(min_ind)
%% Exact solution
% Psi_exact = -1/4/pi./vecnorm(dofs,2,2).^3 .* (dot(dofs,nrms,2));
% 
% % Visualizing the exactness
% figure;
% plot(Psi./Psi_exact);
% 
% err_vec = Psi-Psi_exact;
% 
% % L2 Error
% l2_err = norm(err_vec)
% 
% av_err = err_vec' * V * err_vec
%% Torque evaluation using maxwell stress tensor on the sphere
Xcg = [0 0 0];
centers = mesh_sph.ctr;
Xvec = centers-(centers(:,1)==centers(:,1))*Xcg;
normals_sph = mesh_sph.nrm;
torque_sph_mst = 0.5 * sum( mesh_sph.ndv.*Psi_sph.^2.*cross(Xvec,normals_sph) ,1)

%% Torque evaluation using maxwell stress tensor on the Torus
% Xcg = [0 0 0];
% centers = mesh_tor.ctr;
% Xvec = centers-(centers(:,1)==centers(:,1))*Xcg;
% normals_tor = mesh_tor.nrm;
% torque_tor_mst = 0.5 * sum( mesh_tor.ndv.*Psi_tor.^2.*cross(Xvec,normals_tor) ,1)

%%
% Choice of velocity field for force computation
% Input = N X 3, output = N X 3

% Rotational fields about Xcg
Nux = @(X) (sum(X.*X,2)>R*R).* cross([X(:,1)==X(:,1), 0*X(:,1), 0*X(:,1)],X-Xcg);
Nuy = @(X) (sum(X.*X,2)<R*R).* cross([0*X(:,1), X(:,1)==X(:,1), 0*X(:,1)],X-Xcg);
Nuz = @(X) (sum(X.*X,2)<R*R).* cross([0*X(:,1), 0*X(:,1), X(:,1)==X(:,1)],X-Xcg);

% Visualizing the velocity fields
figure;
plot(mesh);
hold on;
vtcs = S0_Gamma.dof;
vels = Nux(vtcs);
quiver3(vtcs(:,1),vtcs(:,2),vtcs(:,3),vels(:,1),vels(:,2),vels(:,3));
title('Perturbation field');

% Definition of the kernel for T2
kernelx = @(x,y,z) sum(z.*(Nux(x) - Nux(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
kernely = @(x,y,z) sum(z.*(Nuy(x) - Nuy(y)), 2)./(sqrt(sum(z.^2,2)).^3)/ (4*pi);
kernelz = @(x,y,z) sum(z.*(Nuz(x) - Nuz(y)), 2)./(sqrt(sum(z.^2,2)).^3)/ (4*pi);

%t2matx = panel_oriented_assembly(mesh,kernelx,S0_Gamma,S0_Gamma);
%t2maty = panel_oriented_assembly(mesh,kernely,S0_Gamma,S0_Gamma);
%t2matz = panel_oriented_assembly(mesh,kernelz,S0_Gamma,S0_Gamma);

%sum(sum(t2matx))

%forcex = dot(Psi,t2matx * Rho)
%forcey = dot(Psi,t2maty * Rho)
%forcez = dot(Psi,t2matz * Rho)
%forcevals = [forcevals; force];

str1 = "torque_sph_tor_";
fname = append(str1,int2str(N));
%save(fname);
close all;
end