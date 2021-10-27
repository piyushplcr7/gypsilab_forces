% Size of the mesh
N = 100;

% Dimensions of mesh?
L = [1 1 1];

% Creating a mesh
m = mshCube(N,L);

% Defining integration domains?
omega = dom(m,4); % 4 gauss points

% Example of manipulating a mesh
m2 = mshCube(N,2*L);
m2 = translate(m2, [5 5 5]);

% Plotting a mest
plot(m2);

% Getting the boundary of the mesh as a mesh object
mbound = mshBoundary(m);

% domain of integration
Gamma = dom(mbound,3);

% Creating the spaces

% Piecewise linear functions
Vh = fem(m,'P1');

% Ravier thomas 'RWG', Nedelec 'Ned' for other types 


% A way to specify the Dirichlet conditions? Vh had an element called dir?
Vh0 = dirichlet(Vh, mbound);

% Matlab anonymous function? X has 3 columns, function returns the values
% vectorially
G = @(X) (sin(X(:,1)) .* cos(X(:,2)) .* X(:,3));

% Creating a function for grad G?
nablaG = cell(3,1);
nablaG{1} = @(X)(cos(X(:,1).*cos(X(:,2)) .*)s

% Creating integrals 
rhs = integral(Omega,grad(Vh0),gradG);

% Bilinear form
a = integral(Omega, grad(Vh0),grad(Vh0));

%solution
vh = a\rhs;

% Creating a finite element restriction for G
%Gh = G(m.vtx); % does not work because we have to remove boundary points?

% m.vtx contains the boundary that's why we do it this way
%Gh = G(Vh0.dof); % dof for p2 will also contain centers
Gh = G(Vh0.unk); % What is this?

% What is the feval function?



