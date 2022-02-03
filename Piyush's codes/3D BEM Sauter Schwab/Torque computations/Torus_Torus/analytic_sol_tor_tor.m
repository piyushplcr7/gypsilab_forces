close all; clc; clear;
addpath(genpath("../../../"));
format long;
global X;
global W;
load('X3','X');
load('W3','W');

rng(10);
A = rand(3,3);
[Q,R] = qr(A);

% Initializing parameters for the problem

Nvals = 100:50:1700;
sz = size(Nvals,2);

l2errs = zeros(sz,1);
averrs = zeros(sz,1);
linferrs = zeros(sz,1);
Tn_nearest = zeros(sz,1);
Tn_plane = zeros(sz,1);

% Torus radii
r1 = 10;
r2 = 3;

N = 60;
for ii = 1:sz
N = Nvals(ii)

% distance between the centers
dist = 40;

% Mesh for the geometry
mesh_in = mshTorus(N,r1,r2);
mesh_out = mesh_in;

% Rotate the outer torus
mesh_out.vtx = mesh_out.vtx * Q;

% Translate the outer torus
N_vtcs = size(mesh_out.vtx,1);
trans = ones(N_vtcs,1) * [dist 0 0];
mesh_out.vtx = mesh_out.vtx + trans;

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
scatter3(r1+r2,0,0);

%%

% Definition of FEM spaces and integration rule
S1_Gamma = fem(mesh,'P1');
S0_Gamma = fem(mesh,'P0');
S0_Gamma_in = fem(mesh_in,'P0');
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
% Cutoff radius for the BC
%R = 1; 
R = dist/2;
%assert(R < Rad + s);
Vin = 10;
Vout = 20;
% Defining the Dirichlet boundary condition
g = @(x) 1/4/pi./vecnorm(x,2,2); % Point charge inside first torus

% Checking the boundary condition
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
%% Checking the convergence of the Neumann trace at a point
testpt = [r1+r2 0 0];
dofs = S0_Gamma.dof;
Ndofs = size(dofs,1);
distances = vecnorm((dofs - ones(Ndofs,1) * testpt),2,2);
[min_d,min_ind] = min(distances);
Tn_nearest(ii) = Psi(min_ind)
min_d
dofs(min_ind,:)

%% Using the plane method
elts = mesh.elt;
vtcs = mesh.vtx;
Nelts = size(elts,1);

dist_elts = zeros(Nelts,1);

for i = 1:Nelts
    dist_elts(i) = dist_plane(testpt,vtcs(elts(i,1),:),vtcs(elts(i,2),:),vtcs(elts(i,3),:));
end

[mind_plane,minind_plane] = min(dist_elts)

Tn_plane(ii) = Psi(minind_plane)

%% Exact solution
err_vec = Psi;

% L2 Error
l2errs(ii) = err_vec' * M * err_vec %norm(err_vec)

averrs(ii) = err_vec' * V * err_vec

linferrs(ii) = max(abs(err_vec))

close all;

end

save('dirichlet_sol_tor_tor.mat',"averrs","l2errs","Tn_plane","Tn_nearest","Nvals","linferrs");