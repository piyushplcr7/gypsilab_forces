% comparison of proj and projSpecial

addpath(genpath("../../../"));
clear; clc; close all;
format long;

N = 50;
 % Cube size and position
    L = 2*[1 1 1];
    T = [5 5 3];
    
    % Cube domain
    %bndmesh = bndmeshCubeTranslated(N,L,T);
    
    % Spherical domain
%     bndmesh = mshSphere(N,1);
%     bndmesh = bndmesh.translate(T);
    
    mesh = mshCube(N,L);
    mesh = mesh.translate(T);
    %mesh = mesh.sub(1);
    bndmesh = mesh.bnd;
    
    % Mesh size
    %hvals(i) = sqrt(mean(bndmesh.ndv,1));
    
    Gamma = dom(bndmesh,3);
    normals = Gamma.qudNrm;

    % Function to determine the interior of the magnet
    %interior = @(X) (vecnorm(X-T,2,2) < 1);

    % Constant Magnetization function
    M = @(X) ones(size(X,1),1) * [1 0 0];

    [X,W] = Gamma.qud;

    Mvals = M(X);
    Mxn = cross(Mvals,normals,2);

    P1 = fem(bndmesh,'P1');
    DIV0 = P1.nxgrad;

    Mxn_coeffs = projSpecial(Mxn,Gamma,DIV0);
    Mxn_reconstructed = reconstruct(Mxn_coeffs,Gamma,DIV0);

    yoyo = proj(Mxn,Gamma,DIV0);
    yoyor = reconstruct(yoyo,Gamma,DIV0);

    plot(bndmesh);
    hold on;
    quiver3wrapper(X,Mxn,'red');
    quiver3wrapper(X,Mxn_reconstructed,'blue');
    quiver3wrapper(X,yoyor,'green');