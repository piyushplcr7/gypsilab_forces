% New transmission problem script

addpath(genpath("../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 2;
mu0 = 2;
vals = 5:9;
Nvals = size(vals,2);
forces_vol = zeros(Nvals,3);
forces_bem = forces_vol;
torques_vol = forces_vol;
torques_bem = torques_vol;
hvals = vals;

for i = 1:Nvals
    N = 2^vals(i);
    disp(N);
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
    
    Gamma = dom(bndmesh,3);
    normals = Gamma.qudNrm;
    
    %% Source
    N_src = N;
    R0 = 2;
    r0 = .5;
    [J,mesh_src] = get_torus_source(N_src,R0,r0);
    omega_src = dom(mesh_src,3);
    
    %% Solving the transmission problem using scalar potential formulation
    [Tnu,Tdu] = rewrittenSolver(Gamma,omega_src,J,mu,mu0);

    %% Computing the force using MST
    P0 = fem(bndmesh,'P0');
    P1 = fem(bndmesh,'P1');
    % Computing the vector potential Hj at Gamma
    [X_bndmesh,W_bndmesh] = Gamma.qud;
    [Y_src,Wsrc] = omega_src.qud;
    NX = size(X_bndmesh,1);
    NY = size(Y_src,1);
    YY = repmat(Y_src,NX,1);
    XX = repelem(X_bndmesh,NY,1);
    JYY = J(YY(:,1),YY(:,2),YY(:,3));
    gradxG = 1/4/pi * (YY-XX)./vecnorm(XX-YY,2,2).^3;
    integrand = cross(gradxG,JYY,2);
    HJ1 = sum(Wsrc.*reshape(integrand(:,1),NY,NX),1)';
    HJ2 = sum(Wsrc.*reshape(integrand(:,2),NY,NX),1)';
    HJ3 = sum(Wsrc.*reshape(integrand(:,3),NY,NX),1)';

    normals = Gamma.qudNrm;
    HJdotnatGamma = dot([HJ1 HJ2 HJ3],normals,2);
    HJtau = [HJ1 HJ2 HJ3] - HJdotnatGamma.*normals;

    Bn = mu0 * (reconstruct(Tnu,Gamma,P0) + HJdotnatGamma);
    Htau = reconstruct(Tdu,Gamma,grad(P1)) + HJtau;

    jmu = mu0-mu;
    jnu = 1/mu0-1/mu;

    fdensity = 0.5 * (jnu * Bn.^2 - jmu * vecnorm(Htau,2,2).^2).*normals;

    forces_mst(i,:) = sum(W_bndmesh.*fdensity,1)

    %% Computing force using SD BEM
    
    [Vel1,DVel1] = getTransVelDVel([1 0 0]);
    [Vel2,DVel2] = getTransVelDVel([0 1 0]);
    [Vel3,DVel3] = getTransVelDVel([0 0 1]);

    f1 = rewrittenSDBEM(Tdu,Tnu,Gamma,omega_src,mu,mu0,Vel1,DVel1,J);
    f2 = rewrittenSDBEM(Tdu,Tnu,Gamma,omega_src,mu,mu0,Vel2,DVel2,J);
    f3 = rewrittenSDBEM(Tdu,Tnu,Gamma,omega_src,mu,mu0,Vel3,DVel3,J);

    forces_bem(i,:) = [f1 f2 f3]

end
