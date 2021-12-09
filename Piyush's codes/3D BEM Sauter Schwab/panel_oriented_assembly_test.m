% Testing panel oriented assembly code
clear; clc;
% Generating a surface mesh for a sphere using N points
N = 25;
R = 1.5;
mesh = mshSphere(N,R);

% Defining FEM space
S0_Gamma = fem(mesh,'P0');
S1_Gamma = fem(mesh,'P1');

% Defining dom objects for integration in Gypsilab
order = 3;
Gamma = dom(mesh,order);

% Defining the kernel for single layer BIO
% KV = @(x,y,z) 1/norm(z)/4./pi;
KV = @(x,y,z) sqrt(1./ sum(z.^2 ,2) ) /4./pi;

% Getting a matrix using Sauter and Schwab based panel oriented assembly
M = panel_oriented_assembly(mesh,KV,S0_Gamma,S0_Gamma);

% Getting the matrix using Gypsilab implementation
Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); % 0 wave number
V = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,Gxy,S0_Gamma);
V = V + 1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'[1/r]',S0_Gamma);

% Computing the spherical harmonics on mesh 