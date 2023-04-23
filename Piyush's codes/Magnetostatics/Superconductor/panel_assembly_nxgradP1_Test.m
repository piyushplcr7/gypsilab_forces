addpath(genpath("../../../"));
clear; clc; close all;

% Panel assembly double layer test
for N = 50:200:1000
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

    P0 = fem(bndmesh,'P0');
    P1 = fem(bndmesh,'P1');

    Kmat = single_layer(Gamma,P1.nxgrad,P1.nxgrad);

    % Double layer kernel
    kernel = @(x,y,z) 1/(4*pi) * 1./vecnorm(z,2,2); 
    Nelt = bndmesh.nelt;

    [ii,jj] = meshgrid(1:Nelt,1:Nelt);

    Kmat_SS = panel_assembly(bndmesh,kernel,P1.nxgrad,P1.nxgrad,ii(:),jj(:));

    norm(Kmat-Kmat_SS)
end