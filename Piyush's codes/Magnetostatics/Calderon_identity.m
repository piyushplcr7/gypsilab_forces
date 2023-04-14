% Magnetostatics transmission problem with a source

addpath(genpath("../../"));
clear; clc; close all;

mui = 1;
mue = 1;
ivals = 5:12;
err1 = 0*ivals;
err2=err1;
hvals = err1;

nrmerr1 = err1;
nrmerr2 = err1;

Nvals = size(ivals,2);

for i = 1:Nvals
N = 2^ivals(i);
%% SOLUTION DOMAIN
% Cube size and position
L = 2*[1 1 1];
T = [5 5 3];

% Cube domain
%bndmesh = bndmeshCubeTranslated(N,L,T);

% Spherical domain
bndmesh = mshSphere(N,1);
bndmesh = bndmesh.translate(T);

%mesh = mshCube(N,L);
%mesh = mesh.translate(T);
%mesh = mesh.sub(1);
%bndmesh = mesh.bnd;

% Mesh size
hvals(i) = sqrt(mean(bndmesh.ndv,1));

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

%% Galerkin matrices First Equation

% For operator A
Amat = single_layer(Gamma,DIV0,DIV0);
% For operator C
Cmat = double_layer_magnetostatics(Gamma,DIV0,DIV);
% Mass matrix
M_div0_ned = mass_matrix(Gamma,DIV0,NED);

%% galerkin matrices Second Equation
% For operator B, we can simply use C'
Bmat = Cmat';
% For operator N
Nmat = -single_layer(Gamma,DIV.div,DIV.div);

Ndiv0 = size(Amat,1); % Equal to NP1
Nned = size(Nmat,1);

%% Constructing RHS from a current source
[TDA,TNA,mesh_src] = getTracesTorusSource(Gamma);

%[TDA,TNA] = getTraces1(Gamma)
%[TDA,TNA] = getTraces2(Gamma);

% Projecting the Dirichlet trace to NED space
TDA_NED_coeffs = proj(TDA,Gamma,NED);

% Projecting the Neumann trace to DIV space
TNA_DIV_coeffs = proj(TNA,Gamma,DIV);
TNA_DIV0_coeffs = proj(TNA,Gamma,DIV0);

%% Testing the Calderon projector

% Correct version vec1 = Amat * TNA_DIV0_coeffs - (0.5 * M_div0_ned + Cmat) * TDA_NED_coeffs;
vec1 = Amat * TNA_DIV0_coeffs - (0.5 * M_div0_ned + Cmat) * TDA_NED_coeffs;

vec2 = -Nmat * TDA_NED_coeffs + (Cmat'- 0.5 * M_div0_ned') * TNA_DIV0_coeffs;

disp(N)

err1(i) = max(abs(vec1))
%nrmerr1(i) = norm(vec1)
err2(i) = max(abs(vec2))
%nrmerr2(i) = norm(vec2)
%norm(Calderon * traces)

end

loglog(hvals,err1,'-^');
mdl1 = fitlm(log(hvals),log(err1));
m1 = mdl1.Coefficients.Estimate(2);
disp(m1);
c1 = mdl1.Coefficients.Estimate(1);
estimate1 = exp(m1 * log(hvals) + c1);
hold on;

loglog(hvals,estimate1,'-.')

cutoff = 3;
loglog(hvals,err2,'-*');
mdl2 = fitlm(log(hvals(cutoff:end)),log(err2(cutoff:end)));
m2 = mdl2.Coefficients.Estimate(2);
c2 = mdl2.Coefficients.Estimate(1);
disp(m2);
estimate2 = exp(m2 * log(hvals(cutoff:end)) + c2);

loglog(hvals(cutoff:end),estimate2,'--')