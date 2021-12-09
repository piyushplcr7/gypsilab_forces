clear; clc;
% Testing panel oriented assembly code

% Generating a surface mesh for a sphere using N points
N = 50;
R = 1;
mesh = mshSphere(N,R);

% Defining FEM space
S0_Gamma = fem(mesh,'P0');
S1_Gamma = fem(mesh,'P1');

% Defining dom objects for integration in Gypsilab
order = 3;
Gamma = dom(mesh,order);

% 

% Defining the kernel for single layer BIO

KV = @(x,y,z) sqrt(1./ sum(z.^2 ,2) ) /4./pi;

% Getting a matrix using Sauter and Schwab based panel oriented assembly
Vss = panel_oriented_assembly(mesh,KV,S0_Gamma,S0_Gamma);

% Getting the matrix using Gypsilab implementation
Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); % 0 wave number
V = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,Gxy,S0_Gamma);
V = V + 1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'[1/r]',S0_Gamma);

% Defining a harmonic polynomial
g = @(x) x(:,1)+x(:,2)+x(:,3);
l=1;
% Vg = g/(2l+1) cf. thm 8.4 Mclean

L = integral(Gamma,S0_Gamma,g);

sol_gypsi = V\L;

sol_ss = Vss \ L ;

sol_exact = g(mesh.ctr)*(2*l+1);

norm(sol_exact - sol_ss)
norm(sol_gypsi - sol_exact)
