clear; clc;
% Testing panel oriented assembly code

for N = 5:10:1000
% Generating a surface mesh for a sphere using N points
%N = 30;
R = 1.5;
mesh = mshSphere(N,R);

% Defining FEM space
S0_Gamma = fem(mesh,'P0');
S1_Gamma = fem(mesh,'P1');

% Defining dom objects for integration in Gypsilab
order = 3;
Gamma = dom(mesh,order);

% Getting the matrix using Gypsilab implementation
Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); % 0 wave number
V = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,Gxy,S0_Gamma);
V = V + 1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'[1/r]',S0_Gamma);

% Getting the matrix using Gypsilab implementation
% Defining the kernel for DL
GradGn = cell(3,1);
GradGn{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
GradGn{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
GradGn{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);

K = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,GradGn,ntimes(S0_Gamma));
% Regularization
K = K +1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'grady[1/r]',ntimes(S0_Gamma));

% Getting the mass matrix
M = integral(Gamma,S0_Gamma,S0_Gamma);

% Defining the solution
g = @(x) x(:,1)+x(:,2)+x(:,3);

g_N = g(S0_Gamma.dof);

sol = V\((0.5 * M + K)*g_N);

sol_exact = sum(mesh.nrm,2);

err= norm(sol-sol_exact)
end