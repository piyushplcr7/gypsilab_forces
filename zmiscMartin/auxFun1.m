function [p,dp] = auxFun1(curve,t,omega,singVtx,singPow)


t = [curve.I(1);t(:);curve.I(2)];
vtx = [curve.x(t(:)),curve.y(t(:)),0*t(:)];

if ~curve.closed
    elt = [1:(size(vtx,1)-1);2:size(vtx,1)]';
else
    elt = [1:(size(vtx,1));[2:size(vtx,1) 1]]';
end


m = msh(vtx,elt);
[~,B] = m.ABCD;
p = primitiveOnMesh(m,B,omega,singVtx,singPow,10);
p = p(1:end-1);
dp = omega(B(1:end-1,:));


end

