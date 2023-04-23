% Torus source check
clear;
clc;
t = -6:0.5:6;
[X,Y,Z] = meshgrid(t,t,t);

N = size(t,2);

X = reshape(X,[N^3 1]);
Y = reshape(Y,[N^3 1]);
Z = reshape(Z,[N^3 1]);

R0 = 4;
r0 = 2;

interior_torus = @(x,y,z) sqrt( (x-x*R0./(sqrt(x.^2+y.^2))).^2 + (y-y*R0./(sqrt(x.^2+y.^2))).^2 + z.^2 )<r0;
torus_tangent = @(x,y,z) interior_torus(x,y,z).*[-y x z*0]./(sqrt(x.^2+y.^2));

J0 = 1;

J = @(x,y,z) J0*interior_torus(x,y,z).*torus_tangent(x,y,z); 

sizes = ones(N^3,1);
inter = interior_torus(X,Y,Z);
scatter3(X,Y,Z,sizes,inter);

% Vector field check

field = torus_tangent(X,Y,Z);
figure;
quiver3(X,Y,Z,field(:,1),field(:,2),field(:,3));