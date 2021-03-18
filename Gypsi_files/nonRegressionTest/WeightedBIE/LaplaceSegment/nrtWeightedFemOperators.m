%% Boundary element setting

close all;
clear all;
clc;

% Curve creation:
c = openline(-1,1); % Creates a SimpleCurve object

% Refined mesh
m = meshCurve(c,400,'varChange',{@cos,[-pi,0]}); % Creates a msh object

% Singular weight and domain of integration
edges = bnd(m);
X1 = edges.vtx(1,:);
X2 = edges.vtx(2,:);
singVtx = [X1;X2]; % Singularities of 1/omega
singPow = [-1/2;-1/2]; % Power law of the singularities
sing = {singVtx,singPow};
% Weight definition :
[X,Y,Z] = FunR3.XYZ;
omega = sqrt(1 - X^2);
domega{1} = -X./omega;
domega{2} = 0*X;
domega{3} = 0*X;

gss = 5;
Gamma = Wdom(m,gss,1/omega,sing); % Creates a Wdom < dom object
Gamma = Gamma.supplyDw(domega);
% Boundary element space
Vh = P1(m); % Creates a Fe < fem object.


%% Generalized eigenvalue problems

% Finite element matrix in the T^s spaces
Iomega_1 = integral(Gamma,Vh,Vh);
omegaDx2 = integral(Gamma,nxgrad(Vh),omega^2,nxgrad(Vh));

% the generalized eigenvalue problem int omega dxu dxv = lambda int uv/omega
% has solutions given by lambda_n = n^2, and u_n = Tn
[~,D] = eig(full(Iomega_1\omegaDx2));
d = sort(diag(D),'ascend');

disp('Eigenvalues of -(omega dx)^2 :')
disp(' ');
disp(d(1:10))
err = norm(d(1:10)-((0:9).^2)',2)/norm(d(1:10),2);
fprintf('\n Relative error: = %s \n\n',num2str(err));

% Finite element matrix in the U^s spaces
Iomega = integral(Gamma,Vh,omega^2,Vh);

% Weighted Laplacian
dxOmega2 = integral(Gamma,omegaDomega(Vh),omegaDomega(Vh));


% the generalized eigenvalue problem int 1/omega * (omega dx omega) u (omega dx omega) v
% = lambda int omega uv
% has solutions given by lambda_n = (n+1)^2, and u_n = Un
[P,D] = eig(full(Iomega\dxOmega2));
d = sort(diag(D),'ascend');
disp('Eigenvalues of -(dx omega)^2 :')
disp(' ');
disp(d(1:10))
err = norm(d(1:10)-((1:10).^2)',2)/norm(d(1:10),2);
fprintf('\n Relative error = %s \n\n',num2str(err));

%% Galerkin solution and order of convergence.

% 1°) T^s weighted Laplace problem
% Solve ( -(omega dx)^2 + Id) u(x) = f(x)     x in [-1,1]
% With f = [1 + (1-x)^2] sin(x) + x cos(x)
% (exact solution u = sin(x))
% Notice that no boundary conditions are required for this problem.

f = (1 + (1 - X^2))*sin(X) + X*cos(X);
uexact = sin(X);

Ns = [20, 40, 80, 160, 320];
eTs = zeros(size(Ns,2),1);
for i = 1:length(Ns)
    N = Ns(i);
    m = meshCurve(c,N,'varChange',{@cos,[-pi,0]});
    Vh = P1(m);
    Gamma = Wdom(m,gss,1/omega,sing);
    Iomega_1 = integral(Gamma,Vh,Vh);
    omegaDx2 = integral(Gamma,nxgrad(Vh),omega^2,nxgrad(Vh));
    
    rhs = integral(Gamma,Vh,f);
    sol = (Iomega_1 + omegaDx2)\rhs;
    
    e2 = sol'*Iomega_1*sol...
        - 2*sol'*integral(Gamma,Vh,uexact)...
        + integral(Gamma, uexact^2);
    eTs(i) = sqrt(e2);
end



figure;

loglog(Ns,eTs,'-o');
hold on
loglog(Ns,Ns.^(-2),'k--');
legend({'Numerical error','O(1/N^2)'})
xlabel('N')
ylabel('L_{1/\omega}^2 error')
title('T^s Weighted Laplace problem')
eoc = diff(log(eTs'))./diff(log(Ns));
disp('Error:')
disp(eTs);
disp('Estimated order of convergence')
disp(eoc);


% 2°) U^s weighted Laplace problem
% Solve -(dx omega)^2 u(x) = f(x)    x in [-1,1]
% With f = 3X cos(X)  + [1 + (1-X^2)] sin(X)
% (exact solution u = sin(x))
% Notice that no boundary conditions are required for this problem

f = (1 + (1 - X^2))*sin(X) + 3*X*cos(X);
uexact = sin(X);

Ns = [20, 40, 80, 160, 320];
eUs = zeros(size(Ns,2),1);
for i = 1:length(Ns)
    N = Ns(i);
    m = meshCurve(c,N,'varChange',{@cos,[-pi,0]});
    Vh = P1(m);
    Gamma = Wdom(m,gss,1/omega,sing);
    Gamma = Gamma.supplyDw(domega);
    Iomega = integral(Gamma,Vh,(1-X^2),Vh);
    
    dxOmega2 = integral(Gamma,omegaDomega(Vh),omegaDomega(Vh));
    
    rhs = integral(Gamma,Vh,omega^2*f);
    sol =  dxOmega2\rhs;
    
    e2 = sol'*Iomega*sol...
        - 2*sol'*integral(Gamma,Vh,(1-X^2)*uexact)...
        + integral(Gamma, (1-X^2)*(uexact^2));
    eUs(i) = sqrt(e2);
end

figure
loglog(Ns,eUs,'-o');
hold on
loglog(Ns,Ns.^(-2),'k--');
legend({'Numerical error','O(1/N^2)'})
xlabel('N')
ylabel('L_{\omega}^2 error')
title('U^s Weighted Laplace problem')
eoc = diff(log(eUs'))./diff(log(Ns));
disp('Error:')
disp(eUs);
disp('Estimated order of convergence')
disp(eoc);


