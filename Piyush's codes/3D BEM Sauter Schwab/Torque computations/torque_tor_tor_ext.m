%function [] = force_cube_cube_ext(N)

close all; clc; clear;
addpath(genpath("../../../"));
format long;
global X;
global W;
load('X3','X');
load('W3','W');

plotting = false;

% Initializing parameters for the problem

% Torus radii
r1 = 10;
r2 = 3;

N = 60;
for N = 20:100:1500
N
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
rng(10);
A = rand(3,3);
[Q,R] = qr(A);
mesh_in.vtx = mesh_in.vtx * Q;

% Join to create the final mesh
mesh = union(mesh_in,mesh_out);

if plotting
    figure;
    plot(mesh); 
    title('Mesh');
end
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
g = @(X) (sqrt(sum(X.^2,2)) > R)* Vout + (sqrt(sum(X.^2,2)) <= R) * Vin;

% Checking the boundary condition
if plotting
    figure;
    plot(mesh);
    hold on;
    PI_g = g(S0_Gamma.dof);
    plot(mesh,PI_g);
    title('Boundary condition');
    colorbar;
end

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
sPsi = 2*Psi;

if plotting
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
%% Torque evaluation using maxwell stress tensor on the inner torus
Xcg = [0 0 0];
dofs = S0_Gamma_in.dof;
Xvec = dofs-(dofs(:,1)==dofs(:,1))*Xcg;
normals = mesh_in.nrm;
torque_mst = 0.5 * sum( mesh_in.ndv .* Psi_in.^2 .* cross(Xvec,normals) ,1)
force_mst = 0.5 * sum( mesh_in.ndv .* Psi_in.^2 .* normals ,1)

%% Force evaluation for double check

%%
% Choice of velocity field for force computation
% Input = N X 3, output = N X 3

% Velocity fields for the far away sphere
%Nux = @(X) (sum(X.*X,2)>R*R).*[X(:,1)==X(:,1), 0 * X(:,1) , 0 * X(:,1)];
%Nuy = @(X) (sum(X.*X,2)>R*R).*[0*X(:,1), X(:,2)==X(:,2) , 0 * X(:,1)];
%Nuz = @(X) (sum(X.*X,2)>R*R).*[0* X(:,1), 0 * X(:,1) , X(:,3)==X(:,3)];

% Velocity fields for the near sphere
Nux = @(X) (sum(X.*X,2)<R*R).* cross([X(:,1)==X(:,1), 0*X(:,1), 0*X(:,1)],X-Xcg);
%Nunux = @(X) (sum(X.*X,2)<R*R).* [X(:,1)==X(:,1), 0*X(:,1), 0*X(:,1)];
Nuy = @(X) (sum(X.*X,2)<R*R).* cross([0*X(:,1), X(:,1)==X(:,1), 0*X(:,1)],X-Xcg);
Nuz = @(X) (sum(X.*X,2)<R*R).* cross([0*X(:,1), 0*X(:,1), X(:,1)==X(:,1)],X-Xcg);

% Visualizing the velocity fields
if plotting 
    figure;
    plot(mesh);
    hold on;
    vtcs = S0_Gamma.dof;
    vels = Nux(vtcs);
    quiver3(vtcs(:,1),vtcs(:,2),vtcs(:,3),vels(:,1),vels(:,2),vels(:,3));
    title('Perturbation field');
end


% Definition of the kernel for T2
kernelx = @(x,y,z) sum(z.*(Nux(x) - Nux(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
%karnalx = @(x,y,z) sum(z.*(Nunux(x) - Nunux(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
kernely = @(x,y,z) sum(z.*(Nuy(x) - Nuy(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
kernelz = @(x,y,z) sum(z.*(Nuz(x) - Nuz(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);

%t2matx = panel_oriented_assembly(mesh,kernelx,S0_Gamma,S0_Gamma);
%t2fopx = panel_oriented_assembly(mesh,karnalx,S0_Gamma,S0_Gamma);
t2maty = panel_oriented_assembly(mesh,kernely,S0_Gamma,S0_Gamma);
%t2matz = panel_oriented_assembly(mesh,kernelz,S0_Gamma,S0_Gamma);

%sum(sum(t2matx))

%torquex = dot(Psi,t2matx * Rho)
%forcex = dot(Psi,t2fopx * Rho)
torquey = dot(Psi,t2maty * Rho)
%torquez = dot(Psi,t2matz * Rho)

str1 = "Torus_Torus_";
fname = append(str1,int2str(N));
%save(fname);
%exit;

end