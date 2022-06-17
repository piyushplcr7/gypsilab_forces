% Testing panel oriented assembly code

% Generating a surface mesh for a sphere using N points
N = 4;
R = 1.5;
mesh = mshSphere(N,R);

% Defining FEM space
S0_Gamma = fem(mesh,'P0');
S1_Gamma = fem(mesh,'P1');

% Defining dom objects for integration in Gypsilab
order = 3;
Gamma = dom(mesh,order);

% Getting a matrix using Sauter and Schwab based panel oriented assembly
M = panel_oriented_assembly_dl_3d(mesh,S0_Gamma,S0_Gamma);

% Getting the matrix using Gypsilab implementation
% Defining the kernel for DL
GradGn = cell(3,1);
GradGn{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
GradGn{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
GradGn{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);

K = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,GradGn,ntimes(S0_Gamma));
% Regularization
K = K +1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'grady[1/r]',ntimes(S0_Gamma));