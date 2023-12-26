% Equivalent charge simplified script

% Script for constant Magnetization
gpuDevice(2);
addpath(genpath("../../../"));
clear; clc; close all;
format long;

mu = 1;
mu0 = mu;
vals = [1 1/2 1/4 1/7.9 1/15.9];
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
forces_bem = forces_mst;
torques_mst = forces_mst;
torques_bem = forces_mst;
hvals = 0*vals;

for i = 1:Nvals
    N = 2^(i+6); % from 128 to 2048  
    disp(N);
    %% SOLUTION DOMAIN
    % Cube size and position
%     L = [2 0.5 0.5];
%     T = [0 0 0];
% 
%     % Cube domain
%     %bndmesh = bndmeshCubeTranslated(N,L,T);
% 
%     % Spherical domain
%     % bndmesh = mshSphere(N,1);
%     % bndmesh = bndmesh.translate(T);
% 
%     mesh = mshCube(N,L);
%     mesh = mesh.translate(T);
% %     %mesh = mesh.sub(1);
%     bndmesh = mesh.bnd;

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

    % Function to determine the interior of the magnet
    %interior = @(X) (vecnorm(X-T,2,2) < 1);
    
    %% Source
    N_src = N;
    R0 = 2;
    r0 = .5;
    [J_orig,mesh_src] = get_torus_source(N_src,R0,r0);
    omega_src = dom(mesh_src,3);
    
    %% Solving the Magnet problem (SP)
    % Constant Magnetization function
    M = @(X) ones(size(X,1),1) * [1 0 0];
    % Modifying J source
    J = @(x,y,z) 1 * J_orig(x,y,z);

    % Solving for exterior traces
    [psi,g] = solveMagnetProblemSimplifiedSP(Gamma,mu,mu,J,omega_src,M);

    %% Computing the B field for verification
%     plot_field_magnet(TdAM,TnAM_RWG,bndmesh,J,omega_src,mu0,interior);
%     figure;
%     plot_field_magnet(TdAJ,TnAJ_RWG,bndmesh,J,omega_src,mu0,interior);

    %% Computing forces and torques from Equivalent charge model
    
    [X,W] = Gamma.qud;
    Mvals = M(X);
    Mdotn = dot(Mvals,normals,2);

    % Obtaining H values on the magnet boundary
    % {H} = {H.n} + Htan
    HJ = compute_vecpot_curl(J,omega_src,X);

    % Htan = grad_T g
    P1 = fem(bndmesh,'P1');
    Htan = reconstruct(g,Gamma,P1.grad);

    % Normal component is stored in the Neumann trace
    P0 = fem(bndmesh,'P0');
    Hnout = reconstruct(psi,Gamma,P0);
    Hnin = Hnout-Mdotn;
    avgHn = 0.5*(Hnout+Hnin);
    avgH = avgHn.*normals + Htan + HJ;
    
    % Computing the surface integral of (M.n) {H}
    forces_mst(i,:) = mu0 * sum(W.* Mdotn .* avgH,1)

    % Computing force using the simplified expression from BIE based formulation
    forces_bem(i,:) = mu0 * sum(W.* Mdotn .* HJ)


    Xcg = [4 0 0];
    r = X-Xcg;
    torques_mst(i,:) = mu0 * sum(W.* Mdotn .* cross(r,avgH,2),1)

    % Torque computation
    [Velr1,DVelr1] = getRotVelDVel([1 0 0],Xcg);
    [Velr2,DVelr2] = getRotVelDVel([0 1 0],Xcg);
    [Velr3,DVelr3] = getRotVelDVel([0 0 1],Xcg);

    tbem1 = PermanentMagnetShapeDerivativeBIESP(Gamma,g,psi,J,omega_src,Velr1,DVelr1,mu0,M);
    tbem2 = PermanentMagnetShapeDerivativeBIESP(Gamma,g,psi,J,omega_src,Velr2,DVelr2,mu0,M);
    tbem3 = PermanentMagnetShapeDerivativeBIESP(Gamma,g,psi,J,omega_src,Velr3,DVelr3,mu0,M);

    torques_bem(i,:) = [tbem1 tbem2 tbem3]

    save("PM_SP_tetra_1.mat","forces_mst","torques_bem","torques_mst","forces_bem","hvals");
end
