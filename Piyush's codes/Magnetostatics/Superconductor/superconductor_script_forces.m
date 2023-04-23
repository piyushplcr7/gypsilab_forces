% Superconductor script

addpath(genpath("../../../"));
clear; clc; close all;
format long;
mui = 1;
mue = 1;
vals = 5:11;
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
torques_mst = forces_mst;
forces_sd = forces_mst;
torques_sd = forces_mst; 

for i = 1:Nvals
    N = 2^vals(i);
    disp(N)
    %for N=50
    %% SOLUTION DOMAIN
    % Cube size and position
    % L = 2*[1 1 1];
    T = [5 5 3];
    
    % Cube domain
    % bndmesh = bndmeshCubeTranslated(N,L,T);
    
    % Spherical domain
    bndmesh = mshSphere(N,1);
    bndmesh = bndmesh.translate(T);
    
    %mesh = mshCube(N,L);
    %mesh = mesh.translate(T);
    %mesh = mesh.sub(1);
    %bndmesh = mesh.bnd;
    
    % Mesh size
    %hvals(i) = sqrt(mean(bndmesh.ndv,1));

    % Dom object
    Gamma = dom(bndmesh,3);
    
    %% Current source
    N_src = N;
    R0 = 2;
    r0 = .5;
    [J,mesh_src] = get_torus_source(N_src,R0,r0);
    omega_src = dom(mesh_src,3);
    
    %% Solving the problem and obtaining the Neumann trace
    [TnA0,TnA] = solve_superconductor(bndmesh,J,mesh_src);
    
    %% Plotting the computed B field
    
    %plot_field(TnA,bndmesh,J,omega_src);
    
    %% Computing forces
    % Coefficients for zero Dirichlet Trace
    TdA = TnA * 0;
    forces_mst(i,:) = MstForceFromA(TdA,TnA,Gamma)'
    
    % Shape Derivative computation of force
    % Translation fields
    Nux = @(X) ones(size(X,1),1)*[1 0 0];
    Nuy = @(X) ones(size(X,1),1)*[0 1 0];
    Nuz = @(X) ones(size(X,1),1)*[0 0 1];
    
%     sd_e1 = -SuperConductorShapeDerivativeT3(bndmesh,TnA,Nux,omega_src,J);
%     sd_e2 = -SuperConductorShapeDerivativeT3(bndmesh,TnA,Nuy,omega_src,J);
%     sd_e3 = -SuperConductorShapeDerivativeT3(bndmesh,TnA,Nuz,omega_src,J);
% 
%     forces_sd(i,:) = [sd_e1 sd_e2 sd_e3]
    
    %% Computing torques
%     torques_mst(i,:) = MstTorqueFromA(TdA,TnA,Gamma,T)';
% 
%     % Shape Derivative computation of torque
%     % Rotational fields
%     Nuxr = @(X) cross(ones(size(X,1),1)*[1 0 0],X-T);
%     Nuyr = @(X) cross(ones(size(X,1),1)*[0 1 0],X-T);
%     Nuzr = @(X) cross(ones(size(X,1),1)*[0 0 1],X-T);
%     
%     sdt_e1 = -SuperConductorShapeDerivativeT3(bndmesh,TnA,Nuxr,omega_src,J);
%     sdt_e2 = -SuperConductorShapeDerivativeT3(bndmesh,TnA,Nuyr,omega_src,J);
%     sdt_e3 = -SuperConductorShapeDerivativeT3(bndmesh,TnA,Nuzr,omega_src,J);


end