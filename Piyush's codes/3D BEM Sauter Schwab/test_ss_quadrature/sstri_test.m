% Testing of the implemented Sauter Schwab quadrature method.
clear; clc;

% Some test kernels
%Gtest = @(x,y,z) 1; % 0.25
Gtest = @(x,y,z) x(1)*x(2)^2;
%Gtest = @(x,y,z) x(1)*x(2);

% Fundamental solution for Laplace 3D, z := y-x
G1 = @(x,y,z) 1/norm(z)/4./pi;

G2 = @(x,y,z) 1/norm(x-y)/norm(z)/4./pi;

% Double layer kernel in 3D
gradGn1 = @(x,y,z) 1/4./pi * 1/norm(z)^3 * dot(y-x,[0;0;1]);

gradGn2 = @(x,y,z) 1/4./pi * 1/norm(z)^3 * dot(z,[0;0;1]);

gradGn3 = @(x,y,z) 1/4./pi * 1/norm(z)^2 * dot(z/norm(z),[0;0;1]);

% P0 reference shape fuction
hbasis = @(x) 1;
% Parameterization for Triangle with points A,B,C
A = [0;0;0];
B = [1;0;0];
C = [1;1;0];
chi = @(xhat) A + [B-A C-B]*xhat;
% Gram
gram = @(x) 1;

% Sauter Schwab integral
%val = sstri_integrate(G1,hbasis,hbasis,chi,chi,gram,gram,"identical");

%val1 = sstri_integrate(gradGn1,hbasis,hbasis,chi,chi,gram,gram,"identical");

%validentical = sstri_integrate(Gtest,hbasis,hbasis,chi,chi,gram,gram,"identical");
%valce = sstri_integrate(Gtest,hbasis,hbasis,chi,chi,gram,gram,"common_edge");
%valcv = sstri_integrate(Gtest,hbasis,hbasis,chi,chi,gram,gram,"common_vertex");

val = sstri_integrate(Gtest,hbasis,hbasis,chi,chi,gram,gram,"far_away");

%%
% Getting the corresponding integral using Gypsilab for confirmation
clc;
% Getting the mesh first
vtx = [A';B';C';];
elt = [1 2 3];
mesh = msh(vtx,elt);
mesh.ndv

S0_Omega = fem(mesh,'P0');
Omega = dom(mesh,3);
plot(mesh);

% Defining the kernel for SL
Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); % 0 wave number
V = 1/4./pi * integral(Omega,Omega,S0_Omega,Gxy,S0_Omega);
V = V + 1/4./pi * regularize(Omega,Omega,S0_Omega,'[1/r]',S0_Omega);

% Defining the kernel for DL
GradG = cell(3,1);
GradG{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
GradG{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
GradG{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);

K = -1/(4*pi)*integral(Omega,Omega,S0_Omega,GradG,ntimes(S0_Omega));
% Regularization
K = K -1/(4*pi)*regularize(Omega,Omega,S0_Omega,'grady[1/r]',ntimes(S0_Omega));
