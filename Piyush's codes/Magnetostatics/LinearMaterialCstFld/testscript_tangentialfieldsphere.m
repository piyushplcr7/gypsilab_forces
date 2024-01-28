% Testscript LMCF SP

% New transmission problem script
delete(gcp('nocreate'));
addpath(genpath("../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 4;
mu0 = 2;
vals = 8:16;
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
forces_mst_recon = forces_mst;
forces_bem = forces_mst;
torques_mst = forces_mst;
torques_mst_recon = torques_mst;
torques_bem = forces_mst;
hvals = 0*vals;
testsdbem = hvals;
testsdmst = testsdbem;

%rng(32);
H0 = [10 3 1];

for i = 1:Nvals
    
    N = 2^vals(i);
    disp(N);
    %% SOLUTION DOMAIN
    % Cube size and position
    L = 2*[1 1 1];
    T = [1 0.5 2];
    
    % Bounding box
    bndmesh_e = mshSphere(N,9);
    bndmesh_e = bndmesh_e.translate([2 2 2]);

    bndmesh_i = getMeshSphere(N);

    % Mesh size
    hvals(i) = sqrt(mean(bndmesh_i.ndv,1));

    Gamma_i = dom(bndmesh_i,3);
    Gamma_e = dom(bndmesh_e,3);
    normals_i = Gamma_i.qudNrm;

    [X_i,~] = Gamma_i.qud;
    
    [Vel,DVel] = getRotVelDVel([1 0 0],[5 5 3]);
    Vels = Vel(X_i);
    veldotn = dot(Vels,normals_i,2);
    max(abs(veldotn))
%     norm(veldotn)

end
