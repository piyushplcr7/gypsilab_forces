addpath(genpath("../../../../"));
clear; 
clc;
format long;

Nvals = 50:50:4000;
sz = size(Nvals,2);
torques_mst = zeros(sz,3); 
torques_bem = zeros(sz,3);
forces_mst = zeros(sz,3);
forces_bem = zeros(sz,3);
hvals = zeros(sz,1);

for i = 1:sz
    disp(Nvals(i));
    % Get the mesh
    [mesh,mesh_in,mesh_out] = sph_tor_mesh(10,3,5,Nvals(i),10);
    hvals(i) = mean(mesh.ndv,1);
    
    % Solve the floating potential problem on mesh
    [Psi,c] = solve_float_pt_ext(mesh,mesh_in,1e2,3,'gypsi');
    
    S0_Gamma = fem(mesh,'P0');
    Op_in = restriction(S0_Gamma,mesh_in);
    %Psi_in = Op_in * Psi;
    Psi_in = Psi(1:mesh_in.nelt);
    
    % Computing the torques and forces using MST
    [torque_mst,force_mst] = compute_mst_forces(mesh_in,[0,0,0],Psi_in);
    torques_mst(i,:) = torque_mst;
    forces_mst(i,:) = force_mst;
    
    % Computing the torques using BEM formula and parallelization
    %[force_bem,torque_bem] = compute_bem_ft_par(mesh,18,[0,0,0],Psi);
    %torques_bem(i,:) = torque_bem;
    %forces_bem(i,:) = force_bem;

    save('sph_tor_combined.mat','Nvals','forces_mst','torques_bem','torques_mst','forces_bem','hvals');
    
end


