clear all;
close all;

%% msh Mesh 
mesh = msh('sqkite_mesh.msh');
figure;
title('Volume mesh');
plot(mesh);
figure;
bnd_mesh = mesh.bnd;
normals = bnd_mesh.nrm;
centers = bnd_mesh.ctr;
title('Boundary mesh');
plot(bnd_mesh);
hold on;


%% Setting up the discrete system

% Quadrature order
order = 3;

% Creating the integration domain
Omega = dom(mesh, order);
Gamma = dom(bnd_mesh, order);

% FEM space without the dirichlet boundary condition
S1_Gamma = fem(bnd_mesh,'P1');


%% Computing gradu.n DFK

% Computing the Neumann trace
Neu_tr = DirichletDFK(mesh,@g);

% Plotting the Neumann trace, lies in S0(Gamma)
figure;
plot(bnd_mesh);
title('Neumann Trace');
hold on;
pts = S0_Gamma.dof;
scatter3(pts(:,1),pts(:,2),Neu_tr);

%% Computing the volume based expression for force
% Getting quadrature weights and nodes
[QP,QW,elmat] = Omega.qud;
% Getting the unknown to quadrature matrix
uqmat = S0_Omega.uqm(Omega);

% Need to define the function Psi \in P1 and then get coefficients of 
% grad Psi \in P0

% Extracting the elements corresponding to the inner boundary
% element centers are obviously ordered according to the elements
bnd_centers = mesh.bnd.ctr;
dist_bnd_centers = vecnorm(bnd_centers,2,2);
inner_elt_idx = find( dist_bnd_centers < 1.2 );
% mesh.sub(idx) creates a mesh containing elements specified by
% mesh.elt(idx)
inner_bnd_mesh = bnd_mesh.sub(inner_elt_idx);
S1_Gamma_in = fem(inner_bnd_mesh,'P1');

%% Computing the integral 1/2 \int_{\Gamma}(gradu.n)^2 n dS

% Neumann trace lies in S0_Gamma, Neu_tr are the corresponding
% coefficients
Gamma_in = dom(inner_bnd_mesh,1);
S0_Gamma_in = fem(inner_bnd_mesh,'P0');
Mbin_00 = integral(Gamma_in,S0_Gamma_in,ntimes(S0_Gamma_in));
% Restricting the Neumann trace to S0_Gamma_in
S0_Gamma = fem(bnd_mesh,'P0');
% Size S0_Gamma_in_Ndof X S0_Gamma_Ndof
opr = restriction(S0_Gamma,inner_bnd_mesh);
Neu_in = opr * Neu_tr;
F_bnd = zeros(3,1);
F_bnd(1) = 0.5 * Neu_in' * Mbin_00{1} * Neu_in;
F_bnd(2) = 0.5 * Neu_in' * Mbin_00{2} * Neu_in;
F_bnd(3) = 0.5 * Neu_in' * Mbin_00{3} * Neu_in;

%% Plotting the new cut off function along with the mesh
figure;
%plot(mesh);
%hold on;
x = -3:0.08:3;
y=x;
[X,Y] = meshgrid(x,y);
surf(X,Y,cutoff(X,Y));
print('cutoff.eps','-depsc2')


%% Function describing the boundary conditions

function Y = g(X)
    g_o = 0; % Outer boundary condition
    g_i = 1; % Inner boundary condition
    logical_outer = vecnorm(X,2,2)>1.2;
    Y = g_o * logical_outer + g_i * not(logical_outer);
    %Y = X(:,1) + X(:,2);
end
















