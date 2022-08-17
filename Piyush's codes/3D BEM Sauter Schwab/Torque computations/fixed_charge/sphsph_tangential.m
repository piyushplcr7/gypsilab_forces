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

R = 2;
Xcgin = [0 0 0];
Xcgout = [5 0 0];

% Rotational fields
 Nuin = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[1 0 0],X-Xcgin);
 Nuout = @(X) (vecnorm(X,2,2)>=R).* cross(ones(size(X,1),1)*[0 1 0],X-Xcgout);

%Nuin = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[1 0 0].* X,X-Xcgin);
%Nuout = @(X) (vecnorm(X,2,2)>=R).* cross(ones(size(X,1),1)*[0 1 0].*X,X-Xcgout);

% Tangential field
Nu = @(X) (vecnorm(X,2,2)<R).* Nuin(X)+ (vecnorm(X,2,2)>=R).* Nuout(X);

for i = 1:sz
    disp(Nvals(i));
    % Get the mesh
    [mesh,mesh_in,mesh_out] = sph_sph_mesh(1,1,3,Nvals(i));
    hvals(i) = mean(mesh.ndv,1);
    
    % Solve the floating potential problem on mesh
    [Psilol,c] = solve_float_pt_ext(mesh,mesh_in,1e2,3,'gypsi','P0');

    N = size(Psilol,1);
    %Psi = rand(N,1);
    Psi = Psilol;
    
    S0_Gamma = fem(mesh,'P0');
    Op_in = restriction(S0_Gamma,mesh_in);
    %Psi_in = Op_in * Psi;
    Psi_in = Psi(1:mesh_in.nelt);

    % Plotting the velocity field
%     dofs = S0_Gamma.dof;
%     vels = Nu(dofs);
%     plot(mesh);
%     hold on;
%     quiver3(dofs(:,1),dofs(:,2),dofs(:,3), vels(:,1),vels(:,2),vels(:,3));
%     
    % Evaluating the shape derivative formula for tangential fields
    kernel = @(x,y,z) sum(z.*(Nu(x) - Nu(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    t2mat = panel_oriented_assembly(mesh,kernel,S0_Gamma,S0_Gamma);
    val = 0.5 * dot(Psi,t2mat*Psi)

    %save('sph_tor_combined.mat','Nvals','forces_mst','torques_bem','torques_mst','forces_bem','hvals');
    
end


