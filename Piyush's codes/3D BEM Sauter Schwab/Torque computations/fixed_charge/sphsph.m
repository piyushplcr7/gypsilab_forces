addpath(genpath("../../../../"));
clear; 
clc;
format long;

Nvals = 3:12;
Nvals = 2.^Nvals;
sz = size(Nvals,2);
torques_mst = zeros(sz,3); 
torques_bem = zeros(sz,3);
forces_mst = zeros(sz,3);
forces_bem = zeros(sz,1);
hvals = zeros(sz,1);

r1 = 1;
r2 = 2;
d= 3;
R = (r1+r2+d)/2;
Qin = 13.5;

analytical_formula = sphsph_analytic_force(Qin,r1,r2,d)

for i = 1:sz
    disp(Nvals(i));
    % Get the mesh
    [mesh,mesh_in,mesh_out] = sph_sph_mesh(r1,r2,d,Nvals(i));
    hvals(i) = mean(mesh.ndv,1);
    
    % Solve the floating potential problem on mesh
    [Psi,c] = solve_float_pt_ext(mesh,mesh_in,Qin,3,'gypsi','P0');
    
    S0_Gamma = fem(mesh,'P0');
    Op_in = restriction(S0_Gamma,mesh_in);
    %Op_out = restriction(S0_Gamma,mesh_out);
    %Psi_in = Op_in * Psi;
    Psi_in = Psi(1:mesh_in.nelt);
    %Psi_out = Op_out * Psi;
    %outvols = mesh_out.ndv;
    %Qout = sum(outvols.*Psi_out);
    
    % Computing the torques and forces using MST
    [torque_mst,force_mst] = compute_mst_forces(mesh_in,[0,0,0],Psi_in);
    %torques_mst(i,:) = torque_mst;
    forces_mst(i,:) = force_mst
    
    % Computing the torques using BEM formula and parallelization
    Nux = @(X) (vecnorm(X,2,2)<R).* ones(size(X,1),1)*[1 0 0];
    kernelx = @(x,y,z) sum(z.*(Nux(x) - Nux(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    t2mat = panel_oriented_assembly(mesh,kernelx,S0_Gamma,S0_Gamma);

    forces_bem(i) = 0.5 * dot(Psi,t2mat*Psi)

    %coulomb_force = 1/(4*pi) * Qin*Qout/(2*R)^2

    save('sph_sph_analytic.mat','Nvals','forces_mst','forces_bem','hvals','analytical_formula');
    
end
