N = 100;
L = [1 1 1];

m = mshCube(N,L);
% plot(m);
Omega = dom(m,4);
mbound = mshBoundary(m);
Gamma = dom(mbound,3);


Wh = fem(m,'P1');

% Vh = fem(m,'RWG');
% Vh = fem(m,'Ned');

Wh0 = dirichlet(Wh,mbound);

G= @(X)(sin(X(:,1)).*cos(X(:,2)).*X(:,3));

Gh = G(Wh.unk);
Gh0 = G(Wh0.unk);

gradG = cell(3,1);
gradG{1} = @(X)(cos(X(:,1)).*cos(X(:,2)).*X(:,3));
gradG{2} = @(X)(-sin(X(:,1)).*cos(X(:,2)).*X(:,3));
gradG{3} = @(X)(cos(X(:,1)).*cos(X(:,2)));

% Variational formulation:
% int nabla u_h nabla v_h = int nabla v_h . nabla G


% Linear form
rhs = integral(Omega,grad(Wh0),gradG);

% Bilinear form:
a = integral(Omega,grad(Wh0),grad(Wh0));

% Solution:
phih0 = a\rhs;
Phih0 = feval(Wh0,phih0,m);
Uh = Phih0 + Gh;

boundaryValues = feval(Wh,Uh,mbound);

plot(mbound,boundaryValues)
% Now compute the normal trace of uh.

% Now something happens to compute grad u . n 
% this is an element 

Zh = fem(mbound,'P0');

% alpha the coordinates, 
% sum alpha^2 x area of triangles x (normal)




