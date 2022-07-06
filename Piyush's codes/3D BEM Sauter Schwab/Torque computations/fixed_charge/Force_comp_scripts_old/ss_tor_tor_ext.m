close all; clc; clear;
addpath(genpath("../../../../"));
format long;
global X;
global W;
load('X3','X');
load('W3','W');

rng(10);
A = rand(3,3);
[Qrot,R] = qr(A);

% Initializing parameters for the problem

Nvals = 100:100:1700;
sz = size(Nvals,2);

torque_mst = zeros(sz,3);
force_mst = zeros(sz,3);
torque_bem = zeros(sz,3);

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

%% Torque evaluation using maxwell stress tensor on the inner torus
Xcg = [0 0 0];
dofs = S0_Gamma_in.dof;
Xvec = dofs-(dofs(:,1)==dofs(:,1))*Xcg;
normals = mesh_in.nrm;
torque_mst(ii,:) = 0.5 * sum( mesh_in.ndv .* Psi_in.^2 .* cross(Xvec,normals) ,1)
force_mst(ii,:) = 0.5 * sum( mesh_in.ndv .* Psi_in.^2 .* normals ,1)

%%
% Choice of velocity field for force computation
% Input = N X 3, output = N X 3

% Velocity fields for the far away sphere
%Nux = @(X) (sum(X.*X,2)>R*R).*[X(:,1)==X(:,1), 0 * X(:,1) , 0 * X(:,1)];
%Nuy = @(X) (sum(X.*X,2)>R*R).*[0*X(:,1), X(:,2)==X(:,2) , 0 * X(:,1)];
%Nuz = @(X) (sum(X.*X,2)>R*R).*[0* X(:,1), 0 * X(:,1) , X(:,3)==X(:,3)];

% Velocity fields for the near sphere
% Cutoff radius
R = 15;
Nux = @(X) (sum(X.*X,2)<R*R).* cross(X(:,1)*[1 0 0],X-Xcg);
Nuy = @(X) (sum(X.*X,2)<R*R).* cross([0*X(:,1), X(:,1)==X(:,1), 0*X(:,1)],X-Xcg);
Nuz = @(X) (sum(X.*X,2)<R*R).* cross([0*X(:,1), 0*X(:,1), X(:,1)==X(:,1)],X-Xcg);

% Kernels
kernelx = @(x,y,z) sum(z.*(Nux(x) - Nux(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
kernely = @(x,y,z) sum(z.*(Nuy(x) - Nuy(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
kernelz = @(x,y,z) sum(z.*(Nuz(x) - Nuz(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
kernels = cell(3,1);
kernels = {kernelx,kernely,kernelz};
torques = cell(3,1);
t2mats = cell(3,1);
% 
parfor i = 1:3
    disp(i)
    t2mats{i} = panel_oriented_assembly(mesh,kernels{i},S0_Gamma,S0_Gamma);
    torques{i} = dot(Psi,t2mats{i}*Psi);
end

% t2matx = panel_oriented_assembly(mesh,kernelx,S0_Gamma,S0_Gamma);
% torquex = 0.5 * dot(Psi,t2matx*Psi);
% t2maty = panel_oriented_assembly(mesh,kernely,S0_Gamma,S0_Gamma);
% torquey = 0.5 * dot(Psi,t2maty*Psi);
% t2matz = panel_oriented_assembly(mesh,kernelz,S0_Gamma,S0_Gamma);
% torquez = 0.5 * dot(Psi,t2matz*Psi);


% Visualizing the velocity fields
figure;
plot(mesh);
hold on;
vtcs = S0_Gamma.dof;
vels = Nux(vtcs);
quiver3(vtcs(:,1),vtcs(:,2),vtcs(:,3),vels(:,1),vels(:,2),vels(:,3));
title('Perturbation field');

end
save('fp_tor_tor.mat',"torque_mst","force_mst","torque_bem","Nvals");