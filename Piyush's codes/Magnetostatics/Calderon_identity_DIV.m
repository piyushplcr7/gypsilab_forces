% Magnetostatics transmission problem with a source

addpath(genpath("../../"));
clear; clc; close all;

mui = 1;
mue = 1;
for N = 50:200:1000

%% SOLUTION DOMAIN
% Cube size and position
L = 2*[1 1 1];
T = [5 5 3];

% Cube domain
bndmesh = bndmeshCubeTranslated(N,L,T);

% Spherical domain
%bndmesh = mshSphere(N,1);
%bndmesh = bndmesh.translate(T);

%mesh = mshCube(N,L);
%mesh = mesh.translate(T);
%mesh = mesh.sub(1);
%bndmesh = mesh.bnd;

% Mesh size
%hvals(i) = sqrt(mean(bndmesh.ndv,1));

%% BEM
Gamma = dom(bndmesh,3);

% BEM spaces
% curl conforming -> Dirichlet trace
NED = fem(bndmesh,'NED'); 
P1 = fem(bndmesh,'P1');
% Div conforming with div0 constraint -> Neumann trace
DIV0 = nxgrad(P1); 
% Kernel of the surface curl operator
Ker_curl = grad(P1); 
% Div conforming space 
DIV = fem(bndmesh,'RWG');

% Galerkin matrices
% For operator A
Amat = single_layer(Gamma,DIV,DIV);
% For operator C
Cmat = double_layer_magnetostatics(Gamma,DIV,DIV);
% For operator B, we can simply use C'
Bmat = Cmat';
% For operator N
Nmat = -single_layer(Gamma,DIV.div,DIV.div);
% Mass matrix
M_div_ned = mass_matrix(Gamma,DIV,NED);


Ndiv0 = size(Amat,1); % Equal to NP1
Nned = size(Nmat,1);

% Building the interior weak calderon projector
Calderon = [Amat -0.5*M_div_ned-Cmat;
            Cmat'-0.5*M_div_ned' -Nmat];

%% Constructing RHS from a current source
N_src = N;
R0 = 2;
r0 = .5;
[J,mesh_src] = get_torus_source(N_src,R0,r0);
omega_src = dom(mesh_src,3);

% Evaluation points on Gamma
[X,Wt,elt2qud] = Gamma.qud;

% Computing the fields on the points X
A = compute_vecpot(J,omega_src,X);
curlA = compute_vecpot_curl(J,omega_src,X);

% Taking traces of the fields at the points X
normals = Gamma.qudNrm;
TDA = A - dot(A,normals,2).*normals;
TNA = cross(curlA,normals,2);

% Projecting the Dirichlet trace to NED space
TDA_NED_coeffs = proj(TDA,Gamma,NED);

% Projecting the Neumann trace to DIV space
TNA_DIV_coeffs = proj(TNA,Gamma,DIV);

%% Testing the Calderon projector

MDD = mass_matrix(Gamma,DIV,DIV);
MNN = mass_matrix(Gamma,NED,NED);

traces = [TNA_DIV_coeffs; TDA_NED_coeffs];

disp(N)
norm(Calderon * traces)

end