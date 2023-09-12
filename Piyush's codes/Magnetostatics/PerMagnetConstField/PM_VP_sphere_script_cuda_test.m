% Script for constant Magnetization
delete(gcp('nocreate'));
addpath(genpath("../../../"));
clear; clc; close all;
format long;

mu = 4;
mu0 = 2;
vals = 5:11;
Nvals = size(vals,2);
sds = zeros(Nvals,1);
sdsw = zeros(Nvals,1);
forces_mst = zeros(Nvals,3);
forces_new = zeros(Nvals,3);
torques_new = forces_new;
forces_bem = forces_mst;
torques_mst = forces_mst;
torques_bem = forces_mst;
hvals = 0*vals;

for i = 1:Nvals
    N = 2^vals(i);
    disp(N);
    %% SOLUTION DOMAIN
    % Cube size and position
    bndmesh = getMeshSphere(N);
%     bndmesh = getMeshCuboid5(N);
    
    % Mesh size
    hvals(i) = sqrt(mean(bndmesh.ndv,1));
    
    Gamma = dom(bndmesh,3);
    normals = Gamma.qudNrm;

    % Function to determine the interior of the magnet
    %interior = @(X) (vecnorm(X-T,2,2) < 1);
    
    %% Solving the Transmission Problem
    % Constant Magnetic Field function
    B0 = @(X) ones(size(X,1),1) * [1 0 0];

    % These are traces from the exterior
    % Psi lies in the space nxgradP1 and g lies in the space NED
    [TnA,TdA] = solveTransProblemConstField(Gamma,mu,mu0,B0);

    % Projecting TnAM from nxgradP1 to RWG
    P1 = fem(bndmesh,'P1');
    RWG = fem(bndmesh,'RWG');
    Psivals = reconstruct(TnA,Gamma,nxgrad(P1));
    TnA_RWG = proj(Psivals,Gamma,RWG);

    %% Computing the B field for verification
%     plot_field_magnet(TdAM,TnAM_RWG,bndmesh,J,omega_src,mu0,interior);
%     figure;
%     plot_field_magnet(TdAJ,TnAJ_RWG,bndmesh,J,omega_src,mu0,interior);

    %% Computing forces and torques from Equivalent current model

    % Force computation
    NED = fem(bndmesh,'NED'); 
    P1 = fem(bndmesh,'P1');
    P0 = fem(bndmesh,'P0');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0 = nxgrad(P1); 
    RWG = fem(bndmesh,'RWG');

    % Adding normal component of B0
    [X,W] = Gamma.qud;
    B_const_vals = B0(X);
    B_const_n = dot(normals,B_const_vals,2);

    % Bn = curlA.n = curlTg, normal field without the constant field
    Bn = reconstruct(TdA,Gamma,NED.curl);
    Bntotal = Bn + B_const_n;

    % Ht = nx(Hxn) = mu0^-1 nxPsi
    Ht = mu0^(-1) * cross(normals,Psivals,2);
    H_const_t = mu0^(-1) * cross(normals, cross(B_const_vals,normals,2) ,2);
    Httotal = Ht + H_const_t;
%     Httotal = vecnorm(Httotal,2,2);

%     Ht = vecnorm(Ht,2,2);

%     [Vel,DVel] = getCosVelDVel(1,2,1,1);
%     sds(i) = ShapeDervTpVol(Gamma,Bntotal,vecnorm(Httotal,2,2),mu0,mu,Vel)

    forces_mst(i,:) = ForceMstTP(Gamma,Bntotal,vecnorm(Httotal,2,2),mu0,mu)
    
    [Vel1,DVel1] = getTransVelDVel([1 0 0]);
    [Vel2,DVel2] = getTransVelDVel([0 1 0]);
    [Vel3,DVel3] = getTransVelDVel([0 0 1]);

%     sdsw(i) = ForceMstTPConstFiend(Gamma,Bn,Ht,mu0,mu,B0,Vel,DVel,TdA)

    forces_new(i,1) = ForceMstTPConstFiend(Gamma,Bn,Ht,mu0,mu,B0,Vel1,DVel1,TdA);
    forces_new(i,2) = ForceMstTPConstFiend(Gamma,Bn,Ht,mu0,mu,B0,Vel2,DVel2,TdA);
    forces_new(i,3) = ForceMstTPConstFiend(Gamma,Bn,Ht,mu0,mu,B0,Vel3,DVel3,TdA)
% 
    Xcg = [4 0 0];
    torques_mst(i,:) = TorqueMstTP(Gamma,Bntotal,vecnorm(Httotal,2,2),mu0,mu,Xcg)
    [Velr1,DVelr1] = getRotVelDVel([1 0 0],Xcg);
    [Velr2,DVelr2] = getRotVelDVel([0 1 0],Xcg);
    [Velr3,DVelr3] = getRotVelDVel([0 0 1],Xcg);

    torques_new(i,1) = ForceMstTPConstFiend(Gamma,Bn,Ht,mu0,mu,B0,Velr1,DVelr1,TdA);
    torques_new(i,2) = ForceMstTPConstFiend(Gamma,Bn,Ht,mu0,mu,B0,Velr2,DVelr2,TdA);
    torques_new(i,3) = ForceMstTPConstFiend(Gamma,Bn,Ht,mu0,mu,B0,Velr3,DVelr3,TdA)
    
end
