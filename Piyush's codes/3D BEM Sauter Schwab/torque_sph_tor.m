close all; clc; clear;
addpath(genpath("../../"));

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

% Mesh for the geometry
mesh_out = mshSphere(N,Rad);
mesh_in = mshTorus(N,r1,r2);

% Translate the outer sphere
tx = 15;
ty = 0;
tz = 0;
dist = norm([tx ty tz]);

N_vtcs = size(mesh_out.vtx,1);
trans = ones(N_vtcs,1) * [tx ty tz];
mesh_out.vtx = mesh_out.vtx + trans;

% Rotate the inner torus
rng(10);
A = rand(3,3);
[Q,R] = qr(A);
%mesh_in.vtx = mesh_in.vtx * Q;

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
R = 8.5; % Cutoff radius, also used for the velocity field
%g = @(x) (sqrt(sum(x.^2,2)) > R)*(2) + (sqrt(sum(x.^2,2)) <= R)*4;

% Using a known solution
g = @(x) 1/4/pi./sqrt(sum(x.^2,2));

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

%% Testing the obtained Neumann trace

Psi_exact = -1/4/pi./(sqrt(sum(dofs.^2,2)).^3) .* dot(dofs,normals,2);

%% Torque evaluation using maxwell stress tensor on the inner object
Xcg = [0 0 0];
centers = mesh_in.ctr;
Xvec = centers-(centers(:,1)==centers(:,1))*Xcg;
normals_in = mesh_in.nrm;
torque_mst = 0.5 * sum( Psi_in.^2.*cross(Xvec,normals_in) ,1)

%%
% Choice of velocity field for force computation
% Input = N X 3, output = N X 3

% Velocity fields for the far away cube
%Nux = @(X) (sum(X.*X,2)>R*R).*[X(:,1)==X(:,1), 0 * X(:,1) , 0 * X(:,1)];
%Nuy = @(X) (sum(X.*X,2)>R*R).*[0*X(:,1), X(:,2)==X(:,2) , 0 * X(:,1)];
%Nuz = @(X) (sum(X.*X,2)>R*R).*[0* X(:,1), 0 * X(:,1) , X(:,3)==X(:,3)];

% Velocity fields for the near cube
Nux = @(X) (sum(X.*X,2)<R*R).* cross([X(:,1)==X(:,1), 0*X(:,1), 0*X(:,1)],X-Xcg);
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
kernelx = @(x,y,z) sum(z.*(Nux(x) - Nux(y)), 2)./(sqrt(sum(z.^2,2)).^3)/ (4*pi);
kernely = @(x,y,z) sum(z.*(Nuy(x) - Nuy(y)), 2)./(sqrt(sum(z.^2,2)).^3)/ (4*pi);
kernelz = @(x,y,z) sum(z.*(Nuz(x) - Nuz(y)), 2)./(sqrt(sum(z.^2,2)).^3)/ (4*pi);

t2matx = panel_oriented_assembly(mesh,kernelx,S0_Gamma,S0_Gamma);
t2maty = panel_oriented_assembly(mesh,kernely,S0_Gamma,S0_Gamma);
t2matz = panel_oriented_assembly(mesh,kernelz,S0_Gamma,S0_Gamma);

sum(sum(t2matx))

forcex = dot(Psi,t2matx * Rho)
forcey = dot(Psi,t2maty * Rho)
forcez = dot(Psi,t2matz * Rho)
%forcevals = [forcevals; force];

str1 = "torque_sph_tor_";
fname = append(str1,int2str(N));
save(fname);

%end