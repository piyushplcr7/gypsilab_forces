function [alpha,beta,gamma,Nabla] = barycentricCoordinates3(O,A,B,C,X)

N = size(O,1);

a = A-O; 
b = B-O;
c = C-O;
a2 = scal3D(a,a);
b2 = scal3D(b,b);
c2 = scal3D(c,c);
ab = scal3D(a,b);
bc = scal3D(b,c);
ac = scal3D(a,c);


d = a2.*b2.*c2 + 2*ab.*bc.*ac - a2.*bc.^2 - b2.*ac.^2 - c2.*ab.^2;

u = X-O;

v{1} = scal3D(u,a);
v{2} =  scal3D(u,b);
v{3} = scal3D(u,c);


M = cell(3);
M{1,1} = (b2.*c2 - bc.^2)./d;
M{1,2} = (ac.*bc - ab.*c2)./d;
M{1,3} = (ab.*bc - ac.*b2)./d;

M{2,1} = (bc.*ac - ab.*c2)./d;
M{2,2} = (a2.*c2 - ac.^2)./d;
M{2,3} = (ac.*ab - a2.*bc)./d;

M{3,1} = (ab.*bc - b2.*ac)./d;
M{3,2} = (ab.*ac - a2.*bc)./d;
M{3,3} = (a2.*b2 - ab.^2)./d;

r = cell(3,1);
for i = 1:3
    r{i} = zeros(N,1);
    for j = 1:3
        r{i} = r{i} + M{i,j}.*v{j};
    end
end

P = cell(3,3);
P{1,1} = a(:,1); P{1,2} = a(:,2); P{1,3} = a(:,3);
P{2,1} = b(:,1); P{2,2} = b(:,2); P{2,3} = b(:,3);
P{3,1} = c(:,1); P{3,2} = c(:,2); P{3,3} = c(:,3);

Nabla = cellMatrixProd(M,P,@(a,b)(a.*b));

alpha = r{1}; beta = r{2}; gamma = r{3};


end

