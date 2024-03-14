% New transmission problem script

addpath(genpath("../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 100;
mu0 = 1;
vals = [1 1/2 1/4 1/7.9 1/15.9];
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
forces_bem = forces_mst;
torques_mst = forces_mst;
torques_bem = torques_mst;
hvals = vals;

for i = 1:Nvals
    N = 2^(i+6); 
    disp(N);
    %% SOLUTION DOMAIN
    % Cube size and position
    L = [3 1 1];
    T = [2 1 3];
    
    % Cube domain
    %bndmesh = bndmeshCubeTranslated(N,L,T);
    
    % Spherical domain
%     bndmesh = mshSphere(N,1);
%     bndmesh = bndmesh.translate(T);
    
    % mesh = mshCube(N,L);
    % mesh = mesh.translate(T);
    % %mesh = mesh.sub(1);
    % bndmesh = mesh.bnd;
    % bndmesh = meshSymTetra;
    % bndmesh = bndmesh.translate([2 1 3]);
    % bndmesh = bndmesh.refine(vals(i));

    tetra_function_name = sprintf('tetra%d', i);
    tetra_function_handle = str2func(tetra_function_name);
    bndmesh = genMeshFromScript(tetra_function_handle);
    bndmesh = bndmesh.translate([2 1 3]);
    
    % Mesh size
    hvals(i) = sqrt(mean(bndmesh.ndv,1));
    
    Gamma = dom(bndmesh,3);
    normals = Gamma.qudNrm;
    
    %% Source
    N_src = floor((N^2)/350);
    R0 = 2;
    r0 = .5;
    [J,mesh_src] = get_torus_source(N_src,R0,r0);
    omega_src = dom(mesh_src,3);
    
    %% Solving the transmission problem using scalar potential formulation
    % These are traces from the exterior
    % Psi lies in the space nxgradP1 and g lies in the space NED
    [Tnu,Tdu] = solveTpScalPotBIE(bndmesh,mu,mu0,J,omega_src);
    
    %% Computing the MST based force and torque
    [Vel1,DVel1] = getTransVelDVel([1 0 0]);
    [Vel2,DVel2] = getTransVelDVel([0 1 0]);
    [Vel3,DVel3] = getTransVelDVel([0 0 1]);

    f1 = shapeDer_ScalPot_Vol_TP(bndmesh,mu,mu0,Tnu,Tdu,J,omega_src,Vel1,DVel1);
    f2 = shapeDer_ScalPot_Vol_TP(bndmesh,mu,mu0,Tnu,Tdu,J,omega_src,Vel2,DVel2);
    f3 = shapeDer_ScalPot_Vol_TP(bndmesh,mu,mu0,Tnu,Tdu,J,omega_src,Vel3,DVel3);

    forces_mst(i,:) = [f1 f2 f3]

    Xcg = [4 0 0];
    [Velr1,DVelr1] = getRotVelDVel([1 0 0],Xcg);
    [Velr2,DVelr2] = getRotVelDVel([0 1 0],Xcg);
    [Velr3,DVelr3] = getRotVelDVel([0 0 1],Xcg);

    t1 = shapeDer_ScalPot_Vol_TP(bndmesh,mu,mu0,Tnu,Tdu,J,omega_src,Velr1,DVelr1);
    t2 = shapeDer_ScalPot_Vol_TP(bndmesh,mu,mu0,Tnu,Tdu,J,omega_src,Velr2,DVelr2);
    t3 = shapeDer_ScalPot_Vol_TP(bndmesh,mu,mu0,Tnu,Tdu,J,omega_src,Velr3,DVelr3);

    torques_mst(i,:) = [t1 t2 t3]

    % BEM based shape derivative
    fbem1 = shapDervTranPrbScalPotBIE(bndmesh,Tdu,Tnu,J,omega_src,Vel1,DVel1,mu0,mu);
    fbem2 = shapDervTranPrbScalPotBIE(bndmesh,Tdu,Tnu,J,omega_src,Vel2,DVel2,mu0,mu);
    fbem3 = shapDervTranPrbScalPotBIE(bndmesh,Tdu,Tnu,J,omega_src,Vel3,DVel3,mu0,mu);

    forces_bem(i,:) = [fbem1 fbem2 fbem3]

    tbem1 = shapDervTranPrbScalPotBIE(bndmesh,Tdu,Tnu,J,omega_src,Velr1,DVelr1,mu0,mu);
    tbem2 = shapDervTranPrbScalPotBIE(bndmesh,Tdu,Tnu,J,omega_src,Velr2,DVelr2,mu0,mu);
    tbem3 = shapDervTranPrbScalPotBIE(bndmesh,Tdu,Tnu,J,omega_src,Velr3,DVelr3,mu0,mu);

    torques_bem(i,:) = [tbem1 tbem2 tbem3]

    save('TP_SP_tetra_1.mat',"forces_mst","torques_mst","forces_bem","torques_bem","hvals");


end
