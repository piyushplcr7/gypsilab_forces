addpath(genpath("../../../"));
clear; clc; close all;
format long;

vals = 1:9;
Nvals = size(vals,2);

N = 2048;

min_diff = 0*vals;
max_diff = 0*vals;

for i = 1:Nvals
    % N = 2^vals(i);
    % disp(N);
    
%     mesh = mshCube(N,L);
%     mesh = mesh.translate(T);
%     bndmesh = mesh.bnd;

    mesh_src = mshSphere(N,1);

    omega_src = dom(mesh_src,3);

    % Constant surface charge density
    rho = @(X) ones(size(X,1),1);

    % Distance at which potential is computed
    R = 2^vals(i);
    disp(R);

    evalmesh = mshSphere(N,R);
    evalpts = evalmesh.vtx;

    [TdNrho,gradxNrho] = computeNewtonPotentialAndGrad(evalpts,omega_src,rho);

    [X,W] = omega_src.qud;

    rhovals = rho(X);

    Q = sum(W.*rhovals,1);

    exact_val_asymptotic = -Q/4/pi/R^3 * evalpts;

    % max and min difference with exact asymptotic val
    min_diff(i) = min(vecnorm(gradxNrho-exact_val_asymptotic,2,2))
    max_diff(i) = max(vecnorm(gradxNrho-exact_val_asymptotic,2,2))
    

    
end
