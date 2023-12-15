% Superconductor script
gpuDevice(2);
addpath(genpath("../../../"));
clear; clc; close all;
format long;
mui = 1;
mue = 1;
vals = [1 1/2 1/4 1/7.9 1/15.9];
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
torques_mst = forces_mst;
forces_bem = forces_mst;
torques_sd = forces_mst; 
torques_bem = torques_mst;
hvals = vals;

for i = 1:Nvals
    N = 2^(i+6); 
    disp(N)
    %for N=50
    %% SOLUTION DOMAIN
    % Cube size and position
    L = [3 1 1];
    T = [2 1 3];
    
    % Cube domain
    % bndmesh = bndmeshCubeTranslated(N,L,T);
    
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
    
    % plot_field(TnA,bndmesh,J,omega_src);
    
    %% Computing forces
    % Coefficients for zero Dirichlet Trace
    TdA = TnA * 0;
    forces_mst(i,:) = MstForceFromA(TdA,TnA,Gamma)'
    
    % Shape Derivative computation of force
    % Translation fields
    Nux = @(X) ones(size(X,1),1)*[1 0 0];
    Nuy = @(X) ones(size(X,1),1)*[0 1 0];
    Nuz = @(X) ones(size(X,1),1)*[0 0 1];
    
    sd_e1 = -SuperConductorShapeDerivativeT3(bndmesh,TnA,Nux,omega_src,J);
    sd_e2 = -SuperConductorShapeDerivativeT3(bndmesh,TnA,Nuy,omega_src,J);
    sd_e3 = -SuperConductorShapeDerivativeT3(bndmesh,TnA,Nuz,omega_src,J);

    forces_bem(i,:) = [sd_e1 sd_e2 sd_e3]
    
    %% Computing torques
    Xcg = [4 0 0];
    torques_mst(i,:) = MstTorqueFromA(TdA,TnA,Gamma,Xcg)'

    % Shape Derivative computation of torque
    % Getting Rotational Vels and DVels
    [Velxr,DVelxr] = getRotVelDVel([1 0 0],Xcg);
    [Velyr,DVelyr] = getRotVelDVel([0 1 0],Xcg);
    [Velzr,DVelzr] = getRotVelDVel([0 0 1],Xcg);

    rvels = cell(3,1);
    Drvels = cell(3,1);
    rvels{1} = Velxr; rvels{2} = Velyr; rvels{3} = Velzr;
    Drvels{1} = DVelxr; Drvels{2} = DVelyr; Drvels{3} = DVelzr;
    ptorque = zeros(3,1);

    % parfor i = 1:3
    %     ptorque(i) = SuperConductorShapeDerivative(bndmesh,TnA,rvels{i},Drvels{i},omega_src,J)
    % end

    sdt_e1 = SuperConductorShapeDerivative(bndmesh,TnA,Velxr,DVelxr,omega_src,J)
    sdt_e2 = SuperConductorShapeDerivative(bndmesh,TnA,Velyr,DVelyr,omega_src,J)
    sdt_e3 = SuperConductorShapeDerivative(bndmesh,TnA,Velzr,DVelzr,omega_src,J)

    torques_bem(i,:) = [sdt_e1 sdt_e2 sdt_e3]

    save("SC_VP_tetra_1.mat","forces_mst","forces_bem","torques_mst","torques_bem","hvals");

end
