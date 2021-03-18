function [r] = scal3D(X,Y)
% X,Y are Nx3 array containing vector coordinates. 
% X and/or Y may be 1x3.

r = X(:,1).*Y(:,1) + X(:,2).*Y(:,2) + X(:,3).*Y(:,3);

end

