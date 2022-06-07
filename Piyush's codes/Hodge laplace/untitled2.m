addpath(genpath("../../"));
clear; 
clc;
format long;

N = 20;
rad = 5;
mesh = mshSphere(N,rad);
order = 3;
Gamma = dom(mesh,order);

% parameters
stp = mesh.stp;
% Wave number, stable?
k = 1/stp(2);
f = (k*340)/(2*pi);

% FEM spaces
S0_Gamma = fem(mesh,'P0');
S1_Gamma = fem(mesh,'P1');
E0_Gamma = fem(mesh,'NED');

% Kernel
Gxy = @(X,Y)femGreenKernel(X,Y,'[exp(ikr)/r]',k);

% Constructing the block matrix

% Top left
A = integral(Gamma,Gamma,curl(E0_Gamma),Gxy,curl(E0_Gamma));
A = A + regularize(Gamma,Gamma,curl(E0_Gamma),'[1/r]',curl(E0_Gamma));
A = A * (-1./4/pi);

% Top right
B = integral(Gamma,Gamma,nx(E0_Gamma),Gxy,curl(S0_Gamma));
B = B + regularize(Gamma,Gamma,nx(E0_Gamma),'[1/r]',curl(S0_Gamma));
B = B * (-1./4/pi);

% Bottom left
C = B';

% Bottom right identity
D = integral(Gamma,Gamma,nx(S0_Gamma),Gxy,nx(S0_Gamma));
D = D + regularize(Gamma,Gamma,nx(S0_Gamma),'[1/r]',nx(S0_Gamma));
D = D * (-1./4/pi) * k^2;

% Block matrix
T = [A B;
     C D];

% Constructing the RHS





% test space first, trial space second
%size(integral(Gamma,Gamma,S0_Gamma,Gxy,S1_Gamma))

