addpath(genpath("../../../../"));
clear; 
clc;
format long;

Nvals = 100:100:4000;
sz = size(Nvals,2);
forces_mst1 = zeros(sz,3);
forces_gypsi11 = zeros(sz,3);
hvals = zeros(sz,1);

for i = 1:sz
    disp(Nvals(i));
    % Get the mesh
    [mesh,mesh_in,mesh_out] = sph_tor_mesh(10,3,5,Nvals(i),10);
    hvals(i) = mean(mesh.ndv,1);
    
    % Solve the floating potential problem on mesh
    [Psi,c] = solve_float_pt_ext_p1(mesh,mesh_in,1e2,3,'gypsi');
    
    S1_Gamma = fem(mesh,'P1');
    Op_in = restriction(S1_Gamma,mesh_in);
    Psi_in = Op_in * Psi;
    psi_in_space = fem(mesh_in,'P1');
    
    % Computing the torques using BEM formula and parallelization
    forces_gypsi11(i,:) = compute_bem_forces_gypsi(mesh_in,Psi_in,psi_in_space,'P1');

    % Evaluation using MST
    [torque_mst,force_mst] = compute_mst_forces(mesh_in,[0,0,0],Psi_in);
    forces_mst1(i,:) = force_mst;
    
    save('sph_tor_data_gypsi1.mat','Nvals','torques_bem','forces_bem','hvals');
end


