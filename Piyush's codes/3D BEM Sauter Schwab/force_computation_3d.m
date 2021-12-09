% Force computation in 3D
close all; clc; clear;
distances = [];
forcevals = [];

addpath(genpath("../../"));

% Initializing parameters for the problem
N = 30;
R = 2;
offset = 4;
distances = [distances; offset];

% Mesh for the geometry
mesh_in = mshSphere(N,R);
mesh_out = mesh_in;
% Translate the inner sphere
%mesh_in.vtx = mesh_in.vtx + R/4;

% Translate the outer sphere
mesh_out.vtx = mesh_out.vtx + offset;
%mesh_out.translate(offset*[1, 1, 1]);

mesh = union(mesh_in,mesh_out);
%mesh = mesh_in;

figure;
plot(mesh); hold on;

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

% Defining the boundary condition
g = @(X) (sqrt(sum(X.^2,2)) > R*R)*(1) + (sqrt(sum(X.^2,2)) <= R*R)*3;

%g = @(X) 0 * sum(X,2);
%g = @(X) X(:,1)+X(:,2)+X(:,3);
% g = @(X) X(:,1);
%g = @(X) 1./sqrt(sum(X.^2,2));

% Constructing the RHS
% g_N = g(S0_Gamma.dof);
g_N = integral(Gamma,S0_Gamma,g);
% PI_g = M\g_N;

% Checking the boundary condition
%PI_g = g(S0_Gamma.dof);
%plot(mesh,PI_g)
% Solving for state solution
% Interior problem
%Psi = V\((0.5*M + K)*(M\g_N));

gradU{1} = @(X)(-X(:,1)./(X(:,1).^2 + X(:,2).^2 +X(:,3).^2 ).^(3/2));
gradU{2} = @(X)(-X(:,2)./(X(:,1).^2 + X(:,2).^2 +X(:,3).^2 ).^(3/2));
gradU{3} = @(X)(-X(:,3)./(X(:,1).^2 + X(:,2).^2 +X(:,3).^2 ).^(3/2));
Psi_exact = M\integral(Gamma,ntimes(S0_Gamma),gradU);

% Exterior problem
% More exact RHS (Martin's way)
Psi = V\((-0.5 * M + K)* (M\g_N));
% Approximating the RHS
%Psi = V\((-0.5 * M + K)* g_N);

Psi_in = Op_in * Psi;
Psi_out = Op_out * Psi;

% Solving the adjoint problem to get the adjoint solution
% More exact way for RHS (Martin)
Rho = V\(-0.5 * g_N);
%Rho = V\(-0.5 * M * g_N);

% Visualizing the Neumann trace
dofs = S0_Gamma.dof;
normals = mesh.nrm;
sPsi = -2*Psi;
quiver3(dofs(:,1),dofs(:,2),dofs(:,3),normals(:,1).*sPsi,normals(:,2).*sPsi,normals(:,3).*sPsi,0);
%quiver3(dofs(:,1),dofs(:,2),dofs(:,3),normals(:,1),normals(:,2),normals(:,3));

% Checking the solutions
%figure;
%plot(Psi_in);
%hold on;
%plot(Psi_out);

figure;
plot(Psi);
hold on;
% centers = mesh.ctr;
%exact_sol = - sum(centers .* normals,2)./(sqrt(sum(centers.^2,2))).^3;
%plot(exact_solution);
legend(['Psi','Exact solution']);
% plot(sum(normals,2));
plot(Psi_exact);
norm(Psi-Psi_exact)
%%
% Classical formula for force evaluation
normals_in = mesh_in.nrm;
classical_force_in = 0.5 * sum( mesh_in.ndv .* Psi_in.^2.* normals_in, 1)

normals_out = mesh_out.nrm;

classical_force_out = 0.5 * sum( mesh_out.ndv .* Psi_out.^2.* normals_out, 1)

%%
% Choice of velocity field for force computation
% Input = N X 3, output = N X 3
%Nu = @(X) (sum(X.*X,2)<=R*R).*[X(:,1)==X(:,1), 0 * X(:,1) , 0 * X(:,1)];
Nu = @(X) (sum(X.*X,2)>R*R).*[X(:,1)==X(:,1), 0 * X(:,1) , 0 * X(:,1)];
%Nu = @(X) (sum(X.*X,2)>R*R).*[X(:,1)==X(:,1), 0 * X(:,1) , 0 * X(:,1)];
%Nu = @(X) (sum(X.*X,2)<=R*R).*[0 * X(:,1), X(:,1)==X(:,1) , 0 * X(:,1)];
%Nu = @(X) (sum(X.*X,2)<=R*R).*[X(:,1)==X(:,1), 0 * X(:,1) , 0 * X(:,1)];
%Nu = @(X) (sum(X.*X,2)<=R*R).*[X(:,1)==X(:,1), 0 * X(:,1) , 0 * X(:,1)];

% Definition of the kernel for T2
kernel = @(x,y,z) sum(z.*(Nu(x) - Nu(y)), 2)./(sqrt(sum(z.^2,2)).^3)/ (4*pi);

t2mat = panel_oriented_assembly(mesh,kernel,S0_Gamma,S0_Gamma);

sum(sum(t2mat))

force = dot(Psi,t2mat * Rho)
forcevals = [forcevals; force];

exit;