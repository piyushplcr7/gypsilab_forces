% Vecpot evaluation check
addpath(genpath("../../"));
clear;clc;
N_src = 500;
R0 = 2;
r0 = .5;
mesh_src = mshTorus(N_src,R0,r0);
omega_src = dom(mesh_src,3);

% Current
interior_torus = @(x,y,z) sqrt( (x-x*R0./(sqrt(x.^2+y.^2))).^2 + (y-y*R0./(sqrt(x.^2+y.^2))).^2 + z.^2 )<r0;
torus_tangent = @(x,y,z) interior_torus(x,y,z).*[-y x z*0]./(sqrt(x.^2+y.^2));
J0 = 1;
J = @(x,y,z) J0*torus_tangent(x,y,z);

t = -4:0.5:4;
[X,Y,Z] = meshgrid(t,t,t);
sz = size(X);
Ntot = sz(1)* sz(2) * sz(3);

X = reshape(X,[Ntot,1]);
Y = reshape(Y,[Ntot,1]);
Z = reshape(Z,[Ntot,1]);

plot(mesh_src);
hold on;

% Computing the curl of the vector potential, or the magnetic field
eval_X = [X Y Z];
B = compute_vecpot_curl(J,omega_src,eval_X);

quiver3(X,Y,Z,B(:,1),B(:,2),B(:,3));

% Computing on a plane
[X,Z] = meshgrid(t,t);
sz = size(X);
Ntot = sz(1) * sz(2);
X = reshape(X,[Ntot,1]);
Z = reshape(Z,[Ntot,1]);
Y = 0*X;

exterior = ~interior_torus(X,Y,Z);
idx = find(exterior);
X = X(idx);
Y = Y(idx);
Z = Z(idx);

figure;
plot(mesh_src);
hold on;

eval_X = [X Y Z];
B = compute_vecpot_curl(J,omega_src,eval_X);
quiver3(X,Y,Z,B(:,1),B(:,2),B(:,3));

quiver(X,Z,B(:,1),B(:,3));