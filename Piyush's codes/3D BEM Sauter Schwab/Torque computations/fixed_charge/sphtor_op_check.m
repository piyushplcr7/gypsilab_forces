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

for i = 1:sz
    %disp(Nvals(i));
    Nvals(i)
    % Get the mesh
    [mesh,mesh_in,mesh_out] = sph_tor_mesh(10,3,5,Nvals(i),10);
    %hvals(i) = mean(mesh.ndv,1);
    
    % Solve the floating potential problem on mesh
    %[Psi,c] = solve_float_pt_ext(mesh,mesh_in,1e2,3,'gypsi');
    
    S0_Gamma = fem(mesh,'P0');
    Psi = rand(S0_Gamma.ndof,1);
    Op_in = restriction(S0_Gamma,mesh_in);
    Psi_in_op = Op_in * Psi;
    Psi_in = Psi(1:mesh_in.nelt);
    
    norm(Psi_in-Psi_in_op)
    
end

% Result -> the artifical restriction I used works fine!!!


