addpath(genpath("../../../../"));
clear; 
clc;
format long;

Nvals = 3:12;
Nvals = 2.^Nvals;

sz = size(Nvals,2);
sd_full = zeros(sz,1); 
sd_approx = zeros(sz,1);

R = 2;
Xcgin = [0 0 0];
Xcgout = [5 0 0];

% Rotational fields
%Nuin = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[1 0 0],X-Xcgin);
%Nuout = @(X) (vecnorm(X,2,2)>=R).* cross(ones(size(X,1),1)*[0 1 0],X-Xcgout);

Nuin = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[1 0 0].* X,X-Xcgin);
Nuout = @(X) (vecnorm(X,2,2)>=R).* cross(ones(size(X,1),1)*[0 1 0].*X,X-Xcgout);

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
    
    % Evluating near field with Gypsi
    Gamma = dom(mesh,7);
    elt2dof = Gamma.msh.elt;
    
    Ndof = size(Gamma.msh.vtx, 1);
    Nelt = size(elt2dof, 1);

    Adj = sparse((1:Nelt)', elt2dof(:, 1), 1, Nelt, Ndof) + ...
    sparse((1:Nelt)', elt2dof(:, 2), 1, Nelt, Ndof) + ...
    sparse((1:Nelt)', elt2dof(:, 3), 1, Nelt, Ndof);


    [I, J, Case] = find(Adj * Adj');

    corr = ~~(Adj * Adj');

    corr = ~~kron(corr, ones(Gamma.gss));

    % Evaluating the shape derivative formula for tangential fields
    kernel = @(x,y,z) sum(z.*(Nu(x) - Nu(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);

    Gxy = @(x,y) sum((y-x).*(Nu(x) - Nu(y)), 2)./(vecnorm((x-y),2,2).^3)/ (4*pi);

    % Evaluating far field with Gypsi
    A0 = integral(Gamma, Gamma, S0_Gamma, Gxy, S0_Gamma, corr);
    
    t2mat = panel_assembly(mesh,kernel,S0_Gamma,S0_Gamma,I,J);
    t2mat = t2mat + A0;
    sd_approx(i) = 0.5 * dot(Psi,t2mat*Psi)

    t2matfull = panel_oriented_assembly(mesh,kernel,S0_Gamma,S0_Gamma);
    
    sd_full(i) = 0.5 * dot(Psi,t2matfull*Psi)

    save('sph_sph_tangential.mat','Nvals','sd_full','sd_approx','hvals');
    
end


