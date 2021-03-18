function [x,w] = gaussQuadMoments(moments)
% Computes the Gaussian quadrature rule associated with a weight w and 
% endpoints a,b, such that moments(i) = int_{a}^b x^i omega(x) dx
% where i = 0..2n - 2

N = length(moments);
n = fix((N+1)/2);


A = hankel(moments(1:n),moments(n:(2*n-1)));
P = chol(A)^(-1)';
for i = 1:n
    P(i,:) = P(i,:)/P(i,i);
end


% P_i = sum_{k} T_{i,k} x^k
% I = \sum_{k = 0}^{i} \sum_{l = 0}^j T_{i,k} (x^l,x^k) T_{j,l}
% TAT' = I 
% A = (T^{-1}
a = zeros(n,1);
b = zeros(n-1,1);

a(1) = moments(2)/moments(1);
for r = 2:n
    PR = P(r,:);
    XPR = multPoly([0 1],PR);
    XPRPR = multPoly(PR,XPR);
    XPRPR = XPRPR(1:2*r);
    PRPR = multPoly(PR,PR); PRPR = PRPR(1:2*r - 1);
    PR_1 = P(r-1,:);
    PR_1PR_1 = multPoly(PR_1,PR_1); PR_1PR_1 = PR_1PR_1(1:2*(r-1) - 1);
    a(r) = (XPRPR*moments(1:2*r))/(PRPR*moments(1:2*r - 1));
    b(r-1) = (PRPR*moments(1:2*r-1))/(PR_1PR_1*moments(1:2*r-3));
end

J = diag(a) + diag(sqrt(b),-1) + diag(sqrt(b),1);
[P,D] = eig(J);
x = diag(D);
for i = 1:n
    P(:,i) = P(:,i)./sqrt(sum(P(:,i).^2));
end
w = P(1,:)'.^2*moments(1);




end

