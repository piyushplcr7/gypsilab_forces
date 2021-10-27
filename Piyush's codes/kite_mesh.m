clear all;
close all;
N_kite = 50; % Points on the boundary of the kite
theta = linspace(0,2*pi,N_kite+1);
theta = theta(1:N_kite);
theta = theta';
kite_bdry = kite(theta);

% Square boundary
square_bdry = [3 3;
               3 -3;
               -3 -3;
               -3 3];
           
% Combining the points
points = [0 0; kite_bdry; square_bdry];
%points = square_bdry;
DT = delaunayTriangulation(points(:,1),points(:,2));

% Final mesh
elt  = DT.ConnectivityList;
vtx  = [DT.Points,zeros(size(DT.Points,1),1)];
mesh_kite = msh(vtx,elt);
figure;
plot(mesh_kite);

%% Function definitions

function y = kite(x)
  y = [.35 * cos(x) + .1625 * cos(2*x) .35 * sin(x);];
end

function Y = circ_to_kite(X)
    x = X(:,1);
    y = X(:,2);
    [theta,rho] = cart2pol(x,y);
    Y = rho .* kite(theta);
end