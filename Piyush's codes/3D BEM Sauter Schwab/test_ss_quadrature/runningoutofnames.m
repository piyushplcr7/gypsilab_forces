%%
% Getting the corresponding integral using Gypsilab for confirmation
clc;
close all;
clear all;
% Parameterization for Triangle with points A,B,C
A = [0;0;0];
B = [0.6;0;0];
C = [1;1.2;0];
% Getting the mesh first
vtx = [A';B';C';];
elt = [1 2 3];
mesh = msh(vtx,elt);

KV = @(x,y,z) sqrt(1./ sum(z.^2 ,2) ) /4./pi;


mesh=refine(mesh,1);
plot(mesh);


S0_Gamma = fem(mesh,'P0');
% Getting a matrix using Sauter and Schwab based panel oriented assembly
M = panel_oriented_assembly(mesh,KV,S0_Gamma,S0_Gamma);

mesh = refine(mesh,0.025);
figure
plot(mesh);

S0_Omega = fem(mesh,'P0');
Omega = dom(mesh,3);


% Defining the kernel for SL
Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); % 0 wave number
V = 1/4./pi * integral(Omega,Omega,S0_Omega,Gxy,S0_Omega);
V = V + 1/4./pi * regularize(Omega,Omega,S0_Omega,'[1/r]',S0_Omega);

sum(sum(V))-sum(sum(M))