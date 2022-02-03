close all; clc; clear;
addpath(genpath("../../../../"));
format long;
global X;
global W;
load('X3','X');
load('W3','W');
% load('X','X');
% load('W','W');

rng(10);
A = rand(3,3);
[Qrot,R] = qr(A);

% Initializing parameters for the problem

Nvals = 100:50:1700;
sz = size(Nvals,2);

torque_mst = zeros(sz,3);
force_mst = zeros(sz,3);
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

% Translate the outer torus
N_vtcs = size(mesh_out.vtx,1);
trans = ones(N_vtcs,1) * [dist 0 0];
mesh_out.vtx = mesh_out.vtx + trans;

% Rotate the inner torus
mesh_in.vtx = mesh_in.vtx * Qrot;

% Join to create the final mesh
mesh = union(mesh_in,mesh_out);

figure;
plot(mesh); 
title('Mesh');

%%

% Definition of FEM spaces and integration rule
S1_Gamma = fem(mesh,'P1');
S0_Gamma = fem(mesh,'P0');
S0_Gamma_in = fem(mesh_in,'P0');
order = 3;
Gamma = dom(mesh,order);
Gamma_D = dom(mesh_in,order);

Op_in = restriction(S0_Gamma,mesh_in);
Op_out = restriction(S0_Gamma,mesh_out);

% Getting the Single Layer matrix V using Sauter Schwab
kernel = @(x,y,z) 1./vecnorm(z,2,2)/ (4*pi);
V = panel_oriented_assembly(mesh,kernel,S0_Gamma,S0_Gamma);

% Getting the Single Layer matrix V using Gypsilab implementation
Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); % 0 wave number
Vg = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,Gxy,S0_Gamma);
Vg = Vg + 1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'[1/r]',S0_Gamma);

norm(V-Vg)

% l vector for fixed charge formulation
l = integral(Gamma_D, S0_Gamma);

% Block matrix
Q = 10;
Vblock = [V l; l' 0];
temp = Vblock\[0*l; Q];
Psi = temp(1:S0_Gamma.ndof);
c = temp(S0_Gamma.ndof+1);

% Psi in
Psi_in = Op_in * Psi;

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

end