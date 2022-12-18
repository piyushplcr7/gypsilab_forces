% Hypersingular test

addpath(genpath("../../"));
clear; clc;

ivals = 5:12;
Nivals = size(ivals,2);
L2errs = zeros(Nivals,1);
Hhalferrs = L2errs*0;

for i = 1:Nivals
    N = 2^ivals(i);

mesh = mshSphere(N,1.1);
% Integration domain
Gamma = dom(mesh,3);
% BEM spaces
P0 = fem(mesh,'P0');
P1 = fem(mesh,'P1');

% Hypersingular matrix (contains kernel)
W = hypersingular_laplace(Gamma,P1,P1);
% Vector to enforce zero mean for functions
B = integral(Gamma,P1);
% LHS Matrix
A = [W B; B' 0];

% Assembling the RHS for a Neumann trace obtained from a known solution
% u = x+y+z
gradu = [1 1 1];
normals = Gamma.qudNrm;
% Neumann trace at the quadrature points 
TN = sum(gradu.*normals,2);

TN_P0 = proj(TN,Gamma,P0);

% Assembling the RHS
M10 = mass_matrix(Gamma,P1,P0);
% Getting the double layer matrix for the RHS
K = double_layer_laplace(Gamma,P0,P1);

rhs = 0.5*M10*TN_P0 - K' * TN_P0;

rhs = [rhs;0];

% Solving for Dirichlet Trace
sol = A\rhs;
TD_sol = sol(1:end-1);

% Computing the exact Dirichlet trace for comparison
% Quadrature points on Gamma
[X,~,~] = Gamma.qud;
% Dirichlet trace at quadrature points
TD = sum(X,2);
% Projecting to P1
TD_P1 = proj(TD,Gamma,P1);

% Error computation
err_P1 = TD_sol-TD_P1;
disp(N);
err_honehalf = err_P1'*W*err_P1;
M11 = mass_matrix(Gamma,P1,P1);
err_L2 = err_P1'*M11*err_P1;
Hhalferrs(i) = err_honehalf
L2errs(i) = err_L2

end





