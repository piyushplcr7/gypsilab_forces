%%
t = -1:0.01:1;
plot(t.*t,t);
t = 0:0.01:2*pi;
plot(sin(t),cos(2*t));
scatter(sin(t),cos(2*t));
scatter(sin(t),cos(2*t)+cos(3*t));
hold on;
scatter(0.5*sin(t),0.5*cos(2*t)+0.5*cos(3*t));

%% Trying some stuff

N = 40;
x = linspace(0,2*pi, N);
points = kite(x);
X = points(1,:);
Y = points(2,:);
plot(X,Y);

%% Using circle to kite mapping
clear all;
close all;
N = 40;
R = 1;
mesh = mshDisk(N,R);
circ_pts = mesh.vtx;
circ_pts = circ_pts(:,1:2);
kite_pts = circ_to_kite(circ_pts);

scatter(circ_pts(:,1),circ_pts(:,2));
hold on;

scatter(kite_pts(:,1),kite_pts(:,2),'red');

%mesh_kite = delaunay(kite_pts);
% Delaunay triangulation
DT = delaunayTriangulation(kite_pts(:,1),kite_pts(:,2));

% Final mesh
elt  = DT.ConnectivityList;
vtx  = [DT.Points,zeros(size(DT.Points,1),1)];
mesh_kite = msh(vtx,elt);
figure;
plot(mesh_kite);

%% Function definitions

function y = kite(x)
  y = [3.5 * cos(x) + 1.625 * cos(2*x) 3.5 * sin(x);];
end

function Y = circ_to_kite(X)
    x = X(:,1);
    y = X(:,2);
    [theta,rho] = cart2pol(x,y);
    Y = rho .* kite(theta);
end