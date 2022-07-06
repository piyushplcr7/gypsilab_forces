addpath(genpath("../../../../"));
clear; 
clc;
format long;

Nvals = 100:100:4000;
sz = size(Nvals,2);
torques_mst = zeros(sz,3); 
torques_bem = zeros(sz,3);
forces_mst = zeros(sz,3);
forces_bem = zeros(sz,3);
hvals = zeros(sz,1);

R = 12;
Xcg = [0,0,0];
% Translation fields
Nux = @(X) (vecnorm(X,2,2)<R).* ones(size(X,1),1)*[1 0 0];
Nuy = @(X) (vecnorm(X,2,2)<R).* ones(size(X,1),1)*[0 1 0];
Nuz = @(X) (vecnorm(X,2,2)<R).* ones(size(X,1),1)*[0 0 1];

% Rotational fields
Nuxr = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[1 0 0],X-Xcg);
Nuyr = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[0 1 0],X-Xcg);
Nuzr = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[0 0 1],X-Xcg);


for i = 1:sz
    disp(Nvals(i));
    % Get the mesh
    [mesh,mesh_in,mesh_out] = cube_tor_mesh(10,3,[10 10 10],Nvals(i),10);
    hvals(i) = mean(mesh.ndv,1);
    
    % Solve the floating potential problem on mesh
    [Psi,c] = solve_float_pt_ext(mesh,mesh_in,1e2,3,'gypsi');
    
    S0_Gamma = fem(mesh,'P0');
    Op_in = restriction(S0_Gamma,mesh_in);
    Psi_in = Op_in * Psi;
    
    % Computing the torques using BEM formula and parallelization
    force_bem_gypsix = compute_bem_forces_gypsi(mesh,Psi,Nux);
    force_bem_gypsiy = compute_bem_forces_gypsi(mesh,Psi,Nuy);
    force_bem_gypsiz = compute_bem_forces_gypsi(mesh,Psi,Nuz);
    
    torque_bem_gypsix = compute_bem_forces_gypsi(mesh,Psi,Nuxr);
    torque_bem_gypsiy = compute_bem_forces_gypsi(mesh,Psi,Nuyr);
    torque_bem_gypsiz = compute_bem_forces_gypsi(mesh,Psi,Nuzr);

    forces_bem(i,:) = [force_bem_gypsix,force_bem_gypsiy,force_bem_gypsiz];
    torques_bem(i,:) = [torque_bem_gypsix,torque_bem_gypsiy,torque_bem_gypsiz];
    
    save('cube_tor_data_gypsi.mat','Nvals','forces_mst','torques_bem','torques_mst','forces_bem','hvals');
end


