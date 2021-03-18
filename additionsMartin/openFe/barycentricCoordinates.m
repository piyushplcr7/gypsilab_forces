function [alpha,beta,Nabla] = barycentricCoordinates(O,A,B,X)

a = A-O; 
b = B-O;
a2 = scal3D(a,a);
b2 = scal3D(b,b);
ab = scal3D(a,b);
d = a2.*b2 - ab.^2;
mat{1,1} = b2./d;
mat{1,2} = -ab./d;
mat{2,1} = mat{1,2};
mat{2,2} = a2./d;


u = X-O;

ua = scal3D(u,a);
ub = scal3D(u,b);

alpha = mat{1,1}.*ua + mat{1,2}.*ub;
beta = mat{2,1}.*ua + mat{2,2}.*ub;


P{1,1} = a(:,1); P{1,2} = a(:,2); P{1,3} = a(:,3);
P{2,1} = b(:,1); P{2,2} = b(:,2); P{2,3} = b(:,3);

Nabla = cellMatrixProd(mat,P,@(a,b)(a.*b));


end

