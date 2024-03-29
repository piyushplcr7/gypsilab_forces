function [m] = WuniformMesh(curve,N,omega,singVtx,singPow,tol)

if ~ exist('tol','var')||isempty(tol)
    tol = 1e-4;
end

Naccu = 10;

I = curve.I;
t = linspace(I(1),I(2),N+1);
vtx = [curve.x(t(:)),curve.y(t(:)),0*t(:)];

if curve.closed
    elt = [1:N+1;[2:N,1]]';
else
    elt = [1:N;2:(N+1)]';
end

m = msh(vtx,elt);
[~,B] = m.ABCD;

p = WprimitiveOnMesh(m,B,omega,singVtx,singPow,Naccu);
S = sum(p);
Ptarget = linspace(0,S,N+1)';
err = norm(cumsum(p(1:end-1)) - Ptarget(2:end-1),2);
nit = 0;
while ( err > tol && nit < 20)
    t = interp1([0;cumsum(p)],t,Ptarget(2:end-1));
    t = [I(1);t;I(2)];
    vtx = [curve.x(t(:)),curve.y(t(:)),0*t(:)];
    m = msh(vtx,elt);
    [~,B] = m.ABCD;
    p = WprimitiveOnMesh(m,B,omega,singVtx,singPow,Naccu);
    err = norm([0;cumsum(p)] - Ptarget,2);
    S = sum(p);
    Ptarget = linspace(0,S,N+1)';
    nit = nit + 1;
end

