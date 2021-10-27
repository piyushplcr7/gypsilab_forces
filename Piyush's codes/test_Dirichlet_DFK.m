% Test for the function Dirichlet_DFK which solves a dirichlet problem
% using the Direct First Kind BIE

clear all;
close all;
% Creating mesh of annular circle
annular_circle1;
vtx = msh.POS;
elt = msh.TRIANGLES(:,1:3);
clear msh;
mesh = msh(vtx,elt);
bnd_mesh = mesh.bnd;
clear mesh;
mesh = bnd_mesh;
clear bnd_mesh;
plot(mesh);

%% Using the function
psi = Dirichlet_DFK(mesh,@g);

%% Function g representing the boundary conditions
function y = g(X)
    g_o = 0; % Outer boundary condition
    g_i = 1; % Inner boundary condition
    logical_outer = vecnorm(X,2,2)>1;
    y = g_o * logical_outer + g_i * not(logical_outer);
end