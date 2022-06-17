% Solving a Dirichlet BVP on a circle

close all;
clear all;

%% Creating a boundary mesh for a 2D circular domain

R = 1.2; % Radius of the circle
% Center assumed to be the origin

% Mesh parameters
N = 50;

% Creating the mesh manually by constructing the vertices, elements
theta = linspace(0, 2*pi, N+1);
theta = theta(2:N+1);

vtx = [ R * cos(theta); R * sin(theta); 0 * theta]'; % Vertices
elt = [1:N; [2:N] 1]'; % Elements
col = zeros(N,1); % Col?

mesh = msh(vtx,elt,col);
plot(mesh);

% This can also be done by making a 2D mesh and then extracting the
% boundary out of it
% meshvol = mshDisk(N,R);
% mesh = meshvol.bnd;

%% Setting up the discrete system

% Quadrature order
order = 3;
% Creating the integration domain
Gamma = dom(mesh, order);

% Specifying the boundary conditions using function g
%g = @(X)sin(X(:,1)) .* cos(X(:,2));
%g = @(X)ones(size(X(:,1)))*0;
g = @(X)sin(X(:,1)-X(:,2)) .* sinh(X(:,1)+X(:,2));

% Defining the FEM spaces
trial_space = fem(mesh,'P1');
test_space = fem(mesh,'P1');
% A different trial space for the boundary condition g?
trial_space_g = fem(mesh,'P1');

% Defining the kernel for SL
Gxy = @(X,Y)femGreenKernel(X,Y,'[log(r)]',0); % 0 wave number
% Defining the kernel for DL
GradGn = cell(3,1);
GradGn{1} = @(X,Y)femGreenKernel(X,Y,'grady[log(r)]1',0);
GradGn{2} = @(X,Y)femGreenKernel(X,Y,'grady[log(r)]2',0);
GradGn{3} = @(X,Y)femGreenKernel(X,Y,'grady[log(r)]3',0);

% LHS
% Evaluating the bilinear form for the Single Layer BIO
V = -1/(2*pi)*integral(Gamma,Gamma,trial_space,Gxy,test_space);
V = V + -1/(2*pi)*regularize(Gamma,Gamma,trial_space,'[log(r)]',test_space);

% RHS
% Evaluating the bilinear form for the Double Layer BIO
K = -1/(2*pi)*integral(Gamma,Gamma,trial_space,GradGn,ntimes(test_space));
% Regularization
K = K -1/(2*pi)*regularize(Gamma,Gamma,trial_space,'grady[log(r)]',ntimes(test_space));

l = integral(Gamma,trial_space,test_space);

% A clean way to get the coefficients g_N such that when combined with the
% basis elements in trial_space_g, we get an approximation of g?
g_N = g(trial_space_g.dof);

% The linear system is V x = (0.5 l + K) g_N;

% Neumann trace solution:
Psi = V \ ((0.5 * l + K) * g_N );

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
l_ex = integral(Gamma,trial_space,g);

K_ex = -1/(2*pi)*integral(Gamma,Gamma,g,GradGn,ntimes(test_space));
% Regularization
K_ex = K_ex -1/(2*pi)*regularize(Gamma,Gamma,g,'grady[log(r)]',ntimes(test_space));

% Neumann trace solution:
Psi_ex = V \ (0.5 * l_ex + K_ex);


%% Solving the adjoint equation

% The LHS uses the same V matrix computed above. Computing the RHS:

RHS_adj = 0.5 ;

%% Martin stuff

% uqm is unknown to quadrature matrix, used on a fem object
% qud is the quadrature points?








