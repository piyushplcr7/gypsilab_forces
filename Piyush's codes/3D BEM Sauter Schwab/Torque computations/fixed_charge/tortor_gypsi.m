addpath(genpath("../../../../"));
clear; 
clc;
format long;

Nvals = 100:100:4000;
sz = size(Nvals,2);
forces_mst0 = zeros(sz,3);
forces_gypsi00 = zeros(sz,3);
forces_gypsi01 = zeros(sz,3);
%forces_mst1 = zeros(sz,3);
forces_gypsi11 = zeros(sz,3);
hvals = zeros(sz,1);
cvals = hvals;
l2norm0 = cvals;
l2norm1 = cvals;

for i = 1:sz
    disp(Nvals(i));
    % Get the mesh
    [mesh,mesh_in,mesh_out] = tor_tor_mesh(10,3,Nvals(i),10);
    hvals(i) = mean(mesh.ndv,1);
    
    % Solve the floating potential problem on mesh
    [Psi0,c0] = solve_float_pt_ext(mesh,mesh_in,1e2,3,'gypsi','P0');
    [Psi1,c1] = solve_float_pt_ext(mesh,mesh_in,1e2,3,'gypsi','P1');
    
    S0_Gamma = fem(mesh,'P0');
    Op_in0 = restriction(S0_Gamma,mesh_in);
    Psi_in0 = Op_in0 * Psi0;
    psi_in_space0 = fem(mesh_in,'P0');
    
    % Computing the forces usinng double layer formula, P0 for solution
    forces_gypsi00(i,:) = compute_bem_forces_gypsi(mesh,mesh_in,Psi0,S0_Gamma,'P0')
    forces_gypsi01(i,:) = compute_bem_forces_gypsi(mesh,mesh_in,Psi0,S0_Gamma,'P1')

    [torque_mst,force_mst] = compute_mst_forces(mesh_in,[0,0,0],Psi_in0);
    forces_mst0(i,:) = force_mst;

    S1_Gamma = fem(mesh,'P1');
    Op_in1 = restriction(S1_Gamma,mesh_in);
    Psi_in1 = Op_in1 * Psi1;
    psi_in_space1 = fem(mesh_in,'P1');
    
    % Computing the forces usinng double layer formula, P1 for solution
    forces_gypsi11(i,:) = compute_bem_forces_gypsi(mesh,mesh_in,Psi1,S1_Gamma,'P1')

    Gamma = dom(mesh,3);

    M00 = integral(Gamma,S0_Gamma,S0_Gamma);
    M11 = integral(Gamma,S1_Gamma,S1_Gamma);

    l2norm0(i) = dot(Psi0,M00*Psi0);
    l2norm1(i) = dot(Psi1,M11*Psi1);
    
    save('tor_tor_data_gypsi.mat','Nvals','forces_mst0','forces_gypsi00','forces_gypsi01','forces_gypsi11','cvals','l2norm0','l2norm1','hvals');
end




