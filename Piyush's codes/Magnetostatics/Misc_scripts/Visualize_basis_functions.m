% Magnetostatics in vector potential formulation

% Solving curl curl A = J in Omega

addpath(genpath("../../"));

clear; clc;

N = 1;
L = [1 1 1];

mesh = mshCube(N,L);

% Getting single tetrahedron
mesh = mesh.sub(1);

plot(mesh);

V = fem(mesh,'NED');

dofs = V.dof
hold on;
scatter3(dofs(:,1),dofs(:,2),dofs(:,3))

omega = dom(mesh,15);

uqmat = V.uqm(omega);

Ndof = V.ndof;

Id = eye(Ndof,Ndof);
uq1 = uqmat{1} * Id;
uq2 = uqmat{2} * Id;
uq3 = uqmat{3} * Id;

% Visualizing the location of quadrature points
[X,W,elt2qud] = omega.qud;

figure;


m = mesh.bnd;
H = patch('Faces',m.elt,'Vertices',m.vtx);
set(H,'FaceVertexCData',0,'FaceColor','flat');
set(H,'FaceAlpha', 0.2);

hold on;
scatter3(X(:,1),X(:,2),X(:,3));

% Visualizing the ith basis function
i = 4;
%vec = [uq1(:,i) uq2(:,i) uq3(:,i)];
quiver3(X(:,1),X(:,2),X(:,3),uq1(:,i),uq2(:,i),uq3(:,i))

% RWG Element (DIV CONFORMING?)

V = fem(mesh,'RWG');

dofs = V.dof
figure;
scatter3(dofs(:,1),dofs(:,2),dofs(:,3))

omega = dom(mesh,15);

uqmat = V.uqm(omega);

Ndof = V.ndof;

Id = eye(Ndof,Ndof);
uq1 = uqmat{1} * Id;
uq2 = uqmat{2} * Id;
uq3 = uqmat{3} * Id;

% Visualizing the location of quadrature points
[X,W,elt2qud] = omega.qud;

figure;


m = mesh.bnd;
H = patch('Faces',m.elt,'Vertices',m.vtx);
set(H,'FaceVertexCData',0,'FaceColor','flat');
set(H,'FaceAlpha', 0.2);

hold on;
scatter3(X(:,1),X(:,2),X(:,3));

% Visualizing the ith basis function
i = 4;
%vec = [uq1(:,i) uq2(:,i) uq3(:,i)];
quiver3(X(:,1),X(:,2),X(:,3),uq1(:,i),uq2(:,i),uq3(:,i))


