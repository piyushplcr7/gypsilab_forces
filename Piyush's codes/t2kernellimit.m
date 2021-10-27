clear all;
close all;
x = -3:.05:3;
y = x;
N = size(x,2);
[X,Y] = meshgrid(x,y);
points = [reshape(X,[N*N,1]),reshape(Y,[N*N,1])];
vals = t2kernel(points);
vals = reshape(vals,[N,N]);
surf(X,Y,vals);
xlabel('X');
ylabel('Y');


function out = t2kernel(X)
% Point corresponds to a row of X
normx2 = vecnorm(X,2,2).^2;
%disp(size(normx2));
%disp(size(X));
N = size(X,1);
out = sum(X.*( nu(X) - nu(zeros(N,2)) ),2)./normx2;
end

function out = nu(X)
% Point corresponds to a row of X
N = size(X,1);
x = X(:,1);
y = X(:,2);
%out = [sin(x).*cos(y) sin(y)];
out = [sin(x).*cos(y) zeros(N,1)];
end

function out = Dnu(X)
% Point corresponds to a row of X
x = X(:,1);
y = X(:,2);
N = size(X,1);
out = zeros(2,2,N);
out(1,1,:) = cos(x).*cos(y);
out(1,2,:) = -sin(x).*sin(y);
out(2,1,:) = zeros(N,1);
out(2,2,:) = zeros(N,1);
end