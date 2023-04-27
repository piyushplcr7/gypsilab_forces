% Ortho matrix test

addpath(genpath("../../../"));
clear; clc; close all;

for N = 50:100:1000
disp(N);
bndmesh = mshSphere(N,1);

Gamma = dom(bndmesh,3);

NED = fem(bndmesh,'NED'); 
P1 = fem(bndmesh,'P1');
% Div conforming with div0 constraint -> Neumann trace
DIV0 = nxgrad(P1); 
% Kernel of the surface curl operator
Ker_curl = grad(P1); 
% Div conforming space 
DIV = fem(bndmesh,'RWG');

ortho_gypsi = single_layer(Gamma,NED,Ker_curl); % Uses my implementation.

% Use SS to compute the matrix
KV = @(x,y,z) 1./vecnorm(z,2,2)/4./pi;
Nelt = bndmesh.nelt;
[ii,jj] = meshgrid(1:Nelt,1:Nelt);

ortho_SS = -panel_assembly(bndmesh,KV,Ker_curl,nx(DIV),ii(:),jj(:));

norm(ortho_gypsi - ortho_SS)

end