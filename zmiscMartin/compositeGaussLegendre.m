function [X,W] = compositeGaussLegendre(N,n,A,B)

u = linspace(A,B,N+1);
a = u(1:end-1)';
b = u(2:end)';

[x,w] = Gauss_Legendre1D(n,0,1);
X = a + (b-a)*x';
W = (b-a)*w';

X = reshape(X',1,n*N)';
W = reshape(W',1,n*N)';




end

