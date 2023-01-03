% Magnetostatics transmission problem with a source

addpath(genpath("../../"));
clear; clc; close all;

mui = 1;
mue = 1;
N = 50;

%% SOLUTION DOMAIN
% Cube size and position
L = 2*[1 1 1];
T = [5 5 3];

% Cube domain
%bndmesh = bndmeshCubeTranslated(N,L,T);

% Spherical domain
bndmesh = mshSphere(N,1);
bndmesh = bndmesh.translate(T);

% Mesh size
%hvals(i) = sqrt(mean(bndmesh.ndv,1));

%% BEM
Gamma = dom(bndmesh,3);

% BEM spaces
NED = fem(bndmesh,'NED'); % curl conforming -> Dirichlet trace
P1 = fem(bndmesh,'P1');
DIV0 = nxgrad(P1); % div conforming, div0 -> Neumann trace
DIV = fem(bndmesh,'RWG');

% Galerkin matrices

% For operator A
Amat = single_layer(Gamma,DIV0,DIV0);
% For operator C
Cmat = double_layer_magnetostatics(Gamma,DIV0,DIV);
% For operator B, we can simply use C'
Bmat = Cmat';
% For operator N
Nmat = -single_layer(Gamma,DIV.div,DIV.div);
% Vector to enforce zero mean for P1 functions
B = integral(Gamma,P1);

Ndiv0 = size(Amat,1);
Nned = size(Nmat,1);

% Block matrix (Need to subtract constants! matrix is singular)
blockmat = [(mue/mui+1)*Amat -2*Cmat B;...
            2*Bmat -(1+mui/mue)*Nmat zeros(Nned,1);...
            B' zeros(1,Nned) 0];

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
TNA_DIV0_coeffs = proj(TNA,Gamma,DIV0);

M_div0_ned = mass_matrix(Gamma,DIV0,NED);
rhs1 =  mue*(M_div0_ned * TDA_NED_coeffs);
rhs2 = mui * (M_div0_ned' * TNA_DIV0_coeffs);

rhs = [rhs1;rhs2;0];

%% Solving the system

sol = blockmat\rhs;

% Extracting the two traces
TN_sol_coeffs = sol(1:Ndiv0);
TD_sol_coeffs = sol(Ndiv0+1:end-1);

%% Checking with the known traces (which are available for mue/mui=1

TN_sol = reconstruct(TN_sol_coeffs,Gamma,DIV0);
TD_sol = reconstruct(TD_sol_coeffs,Gamma,NED);

quiver3wrapper(X,TNA,'blue');
hold on;
quiver3wrapper(X,TN_sol,'red');

figure;
quiver3wrapper(X,TDA,'blue');
hold on;
quiver3wrapper(X,TD_sol,'red');