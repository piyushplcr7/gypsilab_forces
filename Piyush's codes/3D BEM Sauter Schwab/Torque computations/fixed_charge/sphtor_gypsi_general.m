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

R = 18;
Xcg = [0,0,0];

% Translation fields
%Nux = @(X) (vecnorm(X,2,2)<R).* ones(size(X,1),1)*[1 0 0];
Nux = @(X)  {(vecnorm(X,2,2)<R).* ones(size(X,1),1); 
             (vecnorm(X,2,2)<R).* zeros(size(X,1),1); 
             (vecnorm(X,2,2)<R).* zeros(size(X,1),1); }


for i = 1:sz
    disp(Nvals(i));
    % Get the mesh
    [mesh,mesh_in,mesh_out] = sph_tor_mesh(10,3,5,Nvals(i),10);
    hvals(i) = mean(mesh.ndv,1);
    
    % Solve the floating potential problem on mesh
    [Psi,c] = solve_float_pt_ext(mesh,mesh_in,1e2,3,'gypsi');
    
    S0_Gamma = fem(mesh,'P0');
    Op_in = restriction(S0_Gamma,mesh_in);
    Psi_in = Op_in * Psi;
    psi_in_space = fem(mesh_in,'P0');
    
    % Computing the torques using BEM formula and parallelization
    forces_bem_gypsi = compute_bem_forces_gypsi_general(mesh_in,Psi_in,psi_in_space,Nux)

    %save('sph_tor_data_gypsi.mat','Nvals','torques_bem','forces_bem','hvals');
end


