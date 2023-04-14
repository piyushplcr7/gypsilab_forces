% Checking the representation formula

addpath(genpath("../../"));
clear; clc; close all;

% Cube domain
%bndmesh = bndmeshCubeTranslated(50, ... 
%                                2*[1,1,1], ... % Size of the cube
%                                [5,5,3]); % Translation vector
bndmesh = mshSphere(50,2);
bndmesh = bndmesh.translate([5,5,5]);

% dom object
Gamma = dom(bndmesh,3);

% Source current
[J,mesh_src] = get_torus_source(50, ...
                                2, ... % Ring radius
                                0.5); % Ring thickness

% Source dom object
omega_src = dom(mesh_src,3);

% Quadrature points on Gamma
[Y,W] = Gamma.qud;

% Computing the vector field due to source at points X
A = compute_vecpot(J,omega_src,Y);
curlA = compute_vecpot_curl(J,omega_src,Y);

% Computing traces
normals = Gamma.qudNrm;
TnA = dot(A,normals,2);
TDA = A-TnA.*normals;
RA = cross(normals,TDA,2);
TNA = cross(curlA,normals,2);

% Evaluating the representation formula
% Kernels
grady = cell(3,1);
grady{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
grady{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
grady{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);

Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); 

% Test pt
testpt = [5.75,4.57,5.12];
testptsX = repelem(testpt,size(Y,1),1);

% Computing gradyG(testpts,X)
gradG = [grady{1}(testptsX,Y),...
         grady{2}(testptsX,Y),...
         grady{3}(testptsX,Y)];

% Evaluating \int_{\Gamma} gradyG(x,y) X RA(y) dSy
val1 = sum(W.*( cross(gradG,RA,2) ) ... % Quadrature
            ,1)';

% Evaluating -\int_{\Gamma} G(x,y) (nXcurlA)(y) dSy
%           = \int_{\Gamma} G(x,y) (curlAXn)(y) dSy
val2 = sum(W.*( Gxy(testptsX,Y).* TNA) ...
            ,1)';

% Evaluating -\int_{\Gamma} gradyG(x,y) TnA(y) dSy
val3 = -sum(W.* ( gradG.* TnA) ...
            ,1)';

RHval = (val1-val2+val3)/(4*pi)
PPval = (val1+val2+val3)/(4*pi)
refval = compute_vecpot(J,omega_src,testpt)



