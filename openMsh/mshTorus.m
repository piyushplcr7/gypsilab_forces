function[m] = mshTorus(N,r1,r2)

N1 = fix(sqrt(r1/r2*N))+1;
N2 = fix(N/N1)+1;

a = N1/N2;
m = mshSquare(N,[a,1]);


phi = 2*pi/a*m.vtx(:,1);
theta = 2*pi*m.vtx(:,2);
r = r1 + r2*sin(theta);
z = r2*cos(theta);

m.vtx = [r.*cos(phi),r.*sin(phi),z];

m = swap(m); % normal must go out.
m = clean(m);


end