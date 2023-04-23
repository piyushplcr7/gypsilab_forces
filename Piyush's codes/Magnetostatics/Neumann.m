% curl curl Neumann problem

% Magnetostatics transmission problem with a source

addpath(genpath("../../"));
clear; clc; close all;

mui = 1;
mue = 1;
for i = 6:6
N = 2^i;
%% SOLUTION DOMAIN
% Cube size and position
L = 2*[1 1 1];
T = [5 5 3];

% Cube domain
%bndmesh = bndmeshCubeTranslated(N,L,T);

% Spherical domain
bndmesh = mshSphere(N,1);
bndmesh = bndmesh.translate(T);

% mesh = mshCube(N,L);
% mesh = mesh.translate(T);
% mesh = mesh.sub(1);
% bndmesh = mesh.bnd;

% Mesh size
%hvals(i) = sqrt(mean(bndmesh.ndv,1));

%% BEM
Gamma = dom(bndmesh,3);

% BEM spaces
% curl conforming -> Dirichlet trace
NED = fem(bndmesh,'NED'); 
P1 = fem(bndmesh,'P1');
% Div conforming with div0 constraint -> Neumann trace
%DIV0 = nxgrad(P1); 
% Kernel of the surface curl operator
gradP1 = grad(P1); 
% Div conforming space 
DIV = fem(bndmesh,'RWG');

% Galerkin matrices
% For operator C
Cmat = double_layer_magnetostatics(Gamma,DIV,DIV);
% For operator B, we can simply use C'
Bmat = Cmat';
% For operator N
Nmat = -single_layer(Gamma,DIV.div,DIV.div);
% Vector to enforce zero mean for P1 functions
B = integral(Gamma,P1);
% Matrix to enforce the orthogonal complement constraint
ortho = single_layer(Gamma,NED,gradP1); % Uses my implementation.
%Careful!
%ortho = single_layer(Gamma,Ker_curl,DIV.nx);
%ortho = -ortho';

Ndiv0 = P1.ndof; % Equal to NP1
Nned = NED.ndof;

% Block matrix (Need to subtract constants! matrix is singular)
blockmat = [-Nmat ortho zeros(Nned,1); ortho' zeros(Ndiv0,Ndiv0) B; zeros(1,Nned) B' 0];

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
%TNA_DIV0_coeffs = proj(TNA,Gamma,DIV0);

M_ned_div = mass_matrix(Gamma,NED,DIV);

% rhs = [ (0.5*M_ned_div - Bmat) * TNA_DIV_coeffs ; zeros(Ndiv0,1); 0 ];
% RHS taking into account the weird sign flip
rhs = [ (-0.5*M_ned_div + Bmat) * TNA_DIV_coeffs ; zeros(Ndiv0,1); 0 ];
%% Checking the identity

disp(N)
norm((-Cmat'+0.5*M_ned_div)*TNA_DIV_coeffs -Nmat*TDA_NED_coeffs)

%% Solving the system

sol = blockmat\rhs;

% Extracting the Dirichlet trace
TD_sol_coeffs = sol(1:Nned);

%% Checking with the known traces (which are available for mue/mui=1

TD_sol = reconstruct(TD_sol_coeffs,Gamma,NED);

figure;
quiver3wrapper(X,TDA,'blue');
hold on;
quiver3wrapper(X,TD_sol,'red');

end