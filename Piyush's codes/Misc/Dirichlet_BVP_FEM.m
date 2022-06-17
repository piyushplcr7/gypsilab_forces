% Solving a Dirichlet BVP on a circle

close all;
clear all;
%opengl('save','software')
% Recursive way by genpath? LOL
%addpath(genpath('pathtothatfolder')) 
% ! -> allows to read bash
% e.g. !cd ..
%% Creating a volume mesh for a 2D circular domain

R = 1.2; % Radius of the circle centered at origin

% Mesh parameters
N = 4;

% Using the in built function to create the mesh
mesh = mshDisk(N,R);
bnd_mesh = mesh.bnd;
plot(mesh);

%% Setting up the discrete system

% Quadrature order
order = 3;

% Creating the integration domain
Omega = dom(mesh, order);

% Specifying the boundary conditions using function g
%g = @(X)sin(X(:,1)) .* cos(X(:,2));
%g = @(X)ones(size(X(:,1)))*0;
g = @(X)sin(X(:,1)-X(:,2)) .* sinh(X(:,1)+X(:,2));

% FEM space without the dirichlet boundary condition
space = fem(mesh,'P1');
% FEM space with the homogeneous Dirichlet BC
space0 = dirichlet(space, bnd_mesh);

% LHS matrix
A = integral(Omega,grad(space0),grad(space0));

% How to get the RHS vector? assuming f = -Delta G
% v = -integral(Gamma,f,space0); % Instead of this, use integration by
% parts to convert this expression with only first derivatives of G

% Is it necessary to discretize G? using the extension by zero trick!
% Evaluating g at mesh.vtx would require knowledge of the function on the
% interior so we need zeros at the positions where vertices are in the 
% interior of the mesh

% Initializing to a vector of ones and zeros
%G_N = zeros(size(space.unk,1),1);
%
u_ex = g(space.unk);

% Naive approach
%idx = find(abs(vecnorm(space.unk,2,2)-R) < 1e-12);
%P = zeros(size(space.unk,1),1);
%P(idx) = 1;
%G_N = u_ex .* P;

% Martin's restriction operator approach
P = restriction(space,bnd_mesh); % Gives a matrix with size N_bnd X N_space
g_N = g(bnd_mesh.vtx);
G_N = P' * g_N; % Extension by zero

v = -integral(Omega,grad(space0),grad(space))* G_N;

% This is the offset solution u* 
sol = A\v; 
% Sol contains N_space - N_bnd elements, how to get a solution of size N_space?
% Using Martin's approach of "elimination"
Q = elimination(space,bnd_mesh); % Matrix of size N_space X (N_space-N_bnd)

% Solution with full size N_space, extension by zero
sol = Q * sol;

% Obtaining the solution u by undoing the offset
u = sol + G_N;

err = norm(u-u_ex)

%% Computing the exact solution for testing accuracy

x = R * cos(theta);
y = R * sin(theta);

f_x = cos(x-y) .* sinh(x+y) + sin(x-y) .* cosh(x+y);
f_y = -cos(x-y) .* sinh(x+y) + sin(x-y) .* cosh(x+y);

Psi_exact = f_x .* cos(theta) + f_y .* sin(theta);

figure;

plot(theta,Psi);
hold on;
plot(theta,Psi_exact);

% Exact solution matches with computed solution!

%% Computing without discretizing the boundary condition g
% LHS remains the same, computing RHS:
l_ex = integral(Omega,trial_space,g);

K_ex = -1/(2*pi)*integral(Omega,Omega,g,GradGn,ntimes(test_space));
% Regularization
K_ex = K_ex -1/(2*pi)*regularize(Omega,Omega,g,'grady[log(r)]',ntimes(test_space));

% Neumann trace solution:
Psi_ex = V \ (0.5 * l_ex + K_ex);


%% Solving the adjoint equation

% The LHS uses the same V matrix computed above. Computing the RHS:

RHS_adj = 0.5 ;










