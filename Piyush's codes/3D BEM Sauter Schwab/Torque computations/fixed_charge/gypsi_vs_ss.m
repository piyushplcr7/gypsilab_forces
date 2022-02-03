addpath(genpath("../../../../"));
clear; 
clc;
format long;
global X;
global W;
load('X3','X');
load('W3','W');

Nvals = 50:200:2000;
sz = size(Nvals,2);
torques_mst = zeros(sz,3); 
torques_bem = zeros(sz,3);
forces_mst = zeros(sz,3);
forces_bem = zeros(sz,3);

for i = 1:sz
    
    % Get the mesh
    [mesh,mesh_in,mesh_out] = tor_tor_mesh(10,3,Nvals(i),10);
    
    % Solve the floating potential problem on mesh
    [Psi,c] = solve_float_pt_ext(mesh,mesh_in,1,3,'gypsi');

    % Solve the floating potential problem on mesh
    [Psi1,c1] = solve_float_pt_ext(mesh,mesh_in,1,3,'ss');

    err_psi = norm(Psi-Psi1)
    err_c = abs(c-c1)
end

