delete(gcp('nocreate'));
addpath(genpath("../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
epsilon =1;
epsilon0 = 1;
vals = 5:11;
Nvals = size(vals,2);
energy1 = 0*vals;
energy2 = energy1;

for i = 1:Nvals
    N = 2^vals(i);
    disp(N);
    %% SOLUTION DOMAIN
    % Cube size and position
    L = 2*[1 1 1];
%     T = [1 0.5 2];
    % T = [1 1.5 0.5];
    
%     mesh = mshCube(N,L);
%     mesh = mesh.translate(T);
%     bndmesh = mesh.bnd;

    % Two artificial boundaries
    bndmesh1 = mshSphere(N,1);
    bndmesh2 = mshSphere(N,1);
    bndmesh2 = bndmesh2.translate([-3 -5 -1]);


    %% Source
    mesh_src = mshSphere(N,1);
    % mesh_src = mshCube(N,L);
    mesh_src = mesh_src.translate([3 3 3]);
    omega_src = dom(mesh_src,3);
    % Constant surface charge density
    rho = @(X) ones(size(X,1),1);
    
    %% Solving the transmission problem
    [Tnu1,Tdu1,e1] = solveTpDielSourceChargeEnergy(bndmesh1,epsilon,epsilon0,rho,omega_src);
    [Tnu2,Tdu2,e2] = solveTpDielSourceChargeEnergy(bndmesh2,epsilon,epsilon0,rho,omega_src);
    
    %% Computing energy for the two configurations
    energy1(i) = e1
    energy2(i) = e2
    
end
