% Ortho alternative

addpath(genpath("../../"));
clear; clc; close all;

bndmesh = bndmeshCubeTranslated(50,2*[1 1 1],[5 5 3]);

% BEM spaces
NED = fem(bndmesh,'NED'); % curl conforming -> Dirichlet trace
P1 = fem(bndmesh,'P1');
DIV0 = nxgrad(P1); % div conforming, div0 -> Neumann trace
DIV = fem(bndmesh,'RWG');

Gamma = dom(bndmesh,3);

ortho1 = single_layer(Gamma,NED,P1.grad); % Uses my implementation.

% Ortho 2

kernel = cell(3,1);
kernel{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
kernel{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
kernel{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);

ortho2 = integral(Gamma,Gamma,NED,kernel,P1);
ortho2 = ortho2 + regularize(Gamma,Gamma,NED,'grady[1/r]',P1);
ortho2 = ortho2/(4*pi);

