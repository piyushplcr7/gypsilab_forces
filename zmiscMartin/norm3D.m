function [r] = norm3D(X)
% X is a Nx3 array containing vector coordinates. 
% r(i) is the Euclidean norm of (X(i,1),X(i,2),X(i,3)). 

r = sqrt(scal3D(X,X));

end

