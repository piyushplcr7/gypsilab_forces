% Visualizing basis functions
clear;clc;
eta = 0:0.05:1;
zeta = eta;
[X,Y] = meshgrid(eta,zeta);
X = reshape(X,[21*21,1]);
Y = reshape(Y,[21*21,1]);
inside = (1-X-Y)>=0;
idxinside = find(inside);
X = X(idxinside);
Y = Y(idxinside);
%scatter(X,Y);
%hold on;

basis12 = (1-X-Y).*[1 0] + X.*[1 1];

quiver(X,Y,basis12(:,1),basis12(:,2));

figure;
quiver(X,Y,-basis12(:,1),-basis12(:,2));

figure;
quiver(X,Y,-basis12(:,2),basis12(:,1));

% Visualizing edge elements in 3D
eta = 0:0.1:1;
[X,Y,Z] = meshgrid(eta,eta,eta);
sz = size(X);
N = sz(1)*sz(2)*sz(3);
X = reshape(X,[N,1]);
Y = reshape(Y,[N,1]);
Z = reshape(Z,[N,1]);

inside = (1-X-Y-Z>=0);
idxinside = find(inside);
X = X(idxinside);
Y = Y(idxinside);
Z = Z(idxinside);

basis12 = (1-X-Y-Z).*[1 0 0] - X.*[-1 -1 -1];
figure;
quiver3(X,Y,Z,basis12(:,1),basis12(:,2),basis12(:,3));
xlabel('x');
ylabel('y');
zlabel('z');