%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             BEM for the Hodge-Dirac operator in 3D                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% SPHERE
% Parameters
N = 1e2;
% tol = 1e-4;
k=0;

% Mesh
sphere = mshSphere(N,1);
%S2 = dom(sphere,3);

torus = translate(mshTorus(N,2,1), [4.5,0,0]);
m = union(sphere,torus)
Gamma = dom(m,3)
m.plot



% Green kernel function --> G(x,y) = 1/(4*pi*|x-y|) 
Gxy = @(X,Y) femGreenKernel(X,Y,'[1/r]',k);

% Boundary element spaces
S0 = fem(m, 'P0');
S1 = fem(m, 'P1');
E0x = fem(m, 'RWG');

% Galerkin matrix
A = - 1/(4*pi)*integral(Gamma,Gamma,nxgrad(S1),Gxy,E0x);
A = A - 1/(4*pi)*regularize(Gamma,Gamma,nxgrad(S1),'[1/r]',E0x);
B = -1/(4*pi)*integral(Gamma,Gamma,div(E0x),Gxy,S0);
B = B - 1/(4*pi)*regularize(Gamma,Gamma,div(E0x),'[1/r]',S0);

[n,m] = size(A);
[p,q] = size(B);

Znn = zeros(n,n);
Zmm = zeros(m,m);
Zqq = zeros(q,q);
Znq = zeros(n,q);
Zmp = zeros(m,p);
Zqn = zeros(q,n);

LHS = [Znn A Znq;
       A.' Zmm B;
       Zqn B.' Zqq];

K = null(LHS);
[~, dimK] = size(K);
dimK

%% Torus
% Parameters
N = 1e2;
% tol = 1e-4;
k=0;

% Mesh
torus = mshTorus(N,2,1);
T2 = dom(torus,3);

% Green kernel function --> G(x,y) = 1/(4*pi*|x-y|) 
Gxy = @(X,Y) femGreenKernel(X,Y,'[1/r]',k);

% Boundary element spaces
S0 = fem(torus, 'P0');
S1 = fem(torus, 'P1');
E0x = fem(torus, 'RWG');

% Galerkin matrix
A = - 1/(4*pi)*integral(T2,T2,nxgrad(S1),Gxy,E0x);
A = A - 1/(4*pi)*regularize(T2,T2,nxgrad(S1),'[1/r]',E0x);
B = -1/(4*pi)*integral(T2,T2,div(E0x),Gxy,S0);
B = B - 1/(4*pi)*regularize(T2,T2,div(E0x),'[1/r]',S0);

[n,m] = size(A);
[p,q] = size(B);
n
m
p
q

Znn = zeros(n,n);
Zmm = zeros(m,m);
Zqq = zeros(q,q);
Znq = zeros(n,q);
Zmp = zeros(m,p);
Zqn = zeros(q,n);

LHS = [Znn A Znq;
       A.' Zmm B;
       Zqn B.' Zqq];

K = null(LHS);
[~, dimK] = size(K);
dimK


