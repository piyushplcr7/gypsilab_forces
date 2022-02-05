addpath(genpath("../../../../"));
clear; 
clc;
format long;
global X;
global W;
load('X3','X');
load('W3','W');

Nvals = 50:30:1000;
sz = size(Nvals,2);
torques_mst = zeros(sz,3); 
torques_bem = zeros(sz,3);
forces_mst = zeros(sz,3);
forces_bem = zeros(sz,3);
hvals = zeros(sz,1);

for i = 1:sz
    disp(Nvals(i));
    % Get the mesh
    [mesh,mesh_in,mesh_out] = tor_tor_mesh(10,3,Nvals(i),10);
    hvals(i) = mean(mesh.ndv,1);

    % Solve the floating potential problem on mesh
    [Psi,c] = solve_float_pt_ext(mesh,mesh_in,1,3,'gypsi');
    
    S0_Gamma = fem(mesh,'P0');
    Op_in = restriction(S0_Gamma,mesh_in);
    Psi_in = Op_in * Psi;
    
    % Computing the torques and forces using MST
    [torque_mst,force_mst] = compute_mst_forces(mesh_in,[0,0,0],Psi_in);
    torques_mst(i,:) = torque_mst;
    forces_mst(i,:) = force_mst;
    
    % Computing the torques using BEM formula and parallelization
    fts = zeros(6,1);

    parfor it = 1:6
        switch it
            case 1
                force_bemx = compute_bem_forcex(mesh,18,Psi);
                fts(it) = force_bemx;
            case 2
                force_bemy = compute_bem_forcey(mesh,18,Psi);
                fts(it) = force_bemy;
            case 3
                force_bemz = compute_bem_forcez(mesh,18,Psi);
                fts(it) = force_bemz;
            case 4
                torque_bemx = compute_bem_torquex(mesh,18,[0 0 0],Psi);
                fts(it) = torque_bemx;
            case 5
                torque_bemy = compute_bem_torquey(mesh,18,[0 0 0],Psi);
                fts(it) = torque_bemy;
            case 6
                torque_bemz = compute_bem_torquez(mesh,18,[0 0 0],Psi);
                fts(it) = torque_bemz;
        end
    end
    forces_bem(i,:) = fts(1:3,1)';
    torques_bem(i,:) = fts(4:6,1)';
    
end

%save('tor_tor_data.mat','Nvals','forces_mst','torques_bem','torques_mst','forces_bem','hvals');