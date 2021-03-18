function [r] = norm2D(X)
% X is a Nx2 array containing vector coordinates. 
% r(i) is the Euclidean norm of (X(i,1),X(i,2)). 

r = sqrt(X(:,1).^2 + X(:,2).^2);

end

