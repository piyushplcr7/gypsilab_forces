delete(gcp('nocreate'));
addpath(genpath("../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
epsilon =4;
epsilon0 = 2;
vals = 5:11;
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
forces_mst_recon = forces_mst;
forces_bem = forces_mst;
torques_mst = forces_mst;
torques_mst_recon = torques_mst;
torques_bem = forces_mst;
hvals = 0*vals;

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

    bndmesh = mshSphere(N,1);
    % bndmesh = bndmesh_i.translate([2 0 0]);

    % Mesh size
    hvals(i) = sqrt(mean(bndmesh.ndv,1));

    % Visualize with outer mesh translucent
    %customPlot(bndmesh_i,bndmesh_e);

    %% Source
    mesh_src = mshSphere(N,1);
    % mesh_src = mshCube(N,L);
    mesh_src = mesh_src.translate([3 3 3]);
    omega_src = dom(mesh_src,3);
    % Constant surface charge density
    rho = @(X) ones(size(X,1),1);
    
    %% Solving the transmission problem
    [Tnu,Tdu] = solveTpDielSourceCharge(bndmesh,epsilon,epsilon0,rho,omega_src);
    
    %% Computing force and torque using MST
    Gamma = dom(bndmesh,3);
    P0 = fem(bndmesh,'P0');
    P1 = fem(bndmesh,'P1');

    jumpepsilon = epsilon0 - epsilon;
    jumpInvepsilon = 1/epsilon0 - 1/epsilon;

    grad_g_vals = reconstruct(Tdu,Gamma,grad(P1));
    alpha_vals = epsilon0 * reconstruct(Tnu,Gamma,P0);
    [X,W] = Gamma.qud;
    normals = Gamma.qudNrm;

    fdensity = 0.5 * (jumpepsilon * vecnorm(grad_g_vals,2,2).^2 - jumpInvepsilon * alpha_vals.^2 ) .* normals;
    Xcg = [4 0 0];
    r = X-Xcg;

    forces_mst(i,:) = sum(W.*fdensity,1)
    torques_mst(i,:) = sum(W.* cross(r,fdensity,2),1)

    %% Computing force and torque using BEM SD
    [Vel1,DVel1] = getTransVelDVel([1 0 0]);
    [Vel2,DVel2] = getTransVelDVel([0 1 0]);
    [Vel3,DVel3] = getTransVelDVel([0 0 1]);

    f1 = sdBEMTpDielSC_ConstVel(bndmesh,epsilon,epsilon0,Tnu,Tdu,rho,omega_src,Vel1,DVel1);
    f2 = sdBEMTpDielSC_ConstVel(bndmesh,epsilon,epsilon0,Tnu,Tdu,rho,omega_src,Vel2,DVel2);
    f3 = sdBEMTpDielSC_ConstVel(bndmesh,epsilon,epsilon0,Tnu,Tdu,rho,omega_src,Vel3,DVel3);

    forces_bem(i,:)= [f1 f2 f3]

    [Velr1,DVelr1] = getRotVelDVel([1 0 0],Xcg);
    [Velr2,DVelr2] = getRotVelDVel([0 1 0],Xcg);
    [Velr3,DVelr3] = getRotVelDVel([0 0 1],Xcg);

    t1 = sdBEMTpDielSC(bndmesh,epsilon,epsilon0,Tnu,Tdu,rho,omega_src,Velr1,DVelr1);
    t2 = sdBEMTpDielSC(bndmesh,epsilon,epsilon0,Tnu,Tdu,rho,omega_src,Velr2,DVelr2);
    t3 = sdBEMTpDielSC(bndmesh,epsilon,epsilon0,Tnu,Tdu,rho,omega_src,Velr3,DVelr3);

    torques_bem(i,:) = [t1 t2 t3]
    
end
