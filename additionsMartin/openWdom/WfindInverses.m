function [X,x,FX] = WfindInverses(m,omega,singVtx,singPow,s,X0,tol,BS)

Naccu = 20;
AS = [0; BS(1:end-1)];
g = size(X0,1)/length(m);
S = reshape(s',g,length(m))';
S = S + AS;
S = S';
S = S(:);
tau = repeatLines(m.tgt,g);
X = X0;
A = m.ABCD;
A = repeatLines(A,g);
x = scal3D(X - A,tau);
b = m.ndv;
a = [0; b(1:end-1)];
a = cumsum(a);
x = reshape(x',g,length(m))';
x = x + a;
x = x';
x = x(:);
L = sum(b);

a = repeatLines(a,g);

err = norm(F(X0) - s,2);
for nit = 1:5
    x = interp1([0;F(X);BS(end)],[0;x;L],S);
    X = A + (x - a).*tau;
    FX = F(X);
    err = norm(FX - S,2);
end

for nit = 1:5
    X = X - (FX-S)./omega(X).*tau;
    FX = F(X);
    err = norm(FX - S,2);
end


function[s] = F(X)
s = WprimitiveOnMesh(m,X,omega,singVtx,singPow,Naccu);
s = reshape(s',g,length(m))';
s = s + AS;
s = s';
s = s(:);

end





end

