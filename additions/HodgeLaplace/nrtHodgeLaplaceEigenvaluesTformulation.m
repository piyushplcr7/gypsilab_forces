%% Hodge Laplacian : T formulation


clear all;
close all;

N = 1e2;
mesh = bnd(mshCube(N,[1 1 1]));


Gamma = dom(mesh,3);

Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0);
dG{1} = @(X,Y)femGreenKernel(X,Y,'gradx[1/r]1',0);
dG{2} = @(X,Y)femGreenKernel(X,Y,'gradx[1/r]2',0);
dG{3} = @(X,Y)femGreenKernel(X,Y,'gradx[1/r]3',0);

Vh = fem(mesh,'RWG'); % Vectors Rao-Wilton-Galtson <=> RT0
Wh = fem(mesh,'P0');


% Hodge Laplacian : block matrix of bilinear form
A = 1/(4*pi)*(integral(Gamma,Gamma,Vh,Gxy,Vh) + ...
    regularize(Gamma,Gamma,Vh,'[1/r]',Vh));
B = 1/(4*pi)*(integral(Gamma,Gamma,div(Vh),Gxy,Wh)+ ...
    regularize(Gamma,Gamma,div(Vh),'[1/r]',Wh));
BT = transpose(B);
Z = sparse([],[],[],size(B,2),size(B,2));
system = [[A, B]; [BT, Z]];

% Block mass matrix (in H^-1/2 scalar prod)
M11 = 1/(4*pi)*(integral(Gamma,Gamma,Vh,Gxy,Vh) + ...
    regularize(Gamma,Gamma,Vh,'[1/r]',Vh) + ...
    integral(Gamma,Gamma,div(Vh),Gxy,div(Vh)) + ...
    regularize(Gamma,Gamma,div(Vh),'[1/r]',div(Vh))); 

M22 = 1/(4*pi)*(integral(Gamma,Gamma,Wh,Gxy,Wh) + ...
    regularize(Gamma,Gamma,Wh,'[1/r]',Wh));
M = [M11 0*B; 0*BT M22];


% 0 RHS. 
L = zeros(size(system,1),1);

% 0 solution ?
U = system\L;
disp(max(abs(U)));

% Eigenvalues 
[P,D] = eig(full(M\system));
d = diag(D);

figure;
semilogy(sort(abs(d)),'*');

title('Eigenvalues')

% We find the same result as in Claeys and Hiptmair