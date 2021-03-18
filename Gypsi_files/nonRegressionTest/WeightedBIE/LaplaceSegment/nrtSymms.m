%% Resolution of Symm's integral equation

clear all;
close all;
clc;

%% Mesh and Boundary element space


c = openline(-1,1);
m = meshCurve(c,400,'varChange',{@cos,[-pi,0]});
edges = bnd(m);
X1 = edges.vtx(1,:);
X2 = edges.vtx(2,:);
% Weight definition :
[X,Y,Z] = FunR3.XYZ;
omega1 = sqrt((X1(1) - X).^2 + (X1(2) - Y).^2 + (X1(3) - Z).^2);
omega2 = sqrt((X2(1) - X).^2 + (X2(2) - Y).^2 + (X2(3) - Z).^2);
omega = sqrt(1 - X^2);
singVtx = [X1;X2]; % Singularities of 1/omega
singPow = [-1/2;-1/2]; % Power law of the singularities
sing = {singVtx,singPow};

% Integration domain
gss = 5;
tic;
Gamma = Wdom(m,gss,1/omega,sing);
toc;
check = integral(Gamma,0*X +1);
pi - check
% Finite element:
Vh = P1(m);


%% Generalized eigenvalue problem

k = 0;
GXY = @(X,Y)femGreenKernel(X,Y,'[log(r)]',k);
Somega = -1/(2*pi)*(...
    integral(Gamma,Gamma,Vh,GXY,Vh)  ...
    + regularize(Gamma,Gamma,Vh,'[log(r)]',Vh));

Iomega_1 = integral(Gamma,Vh,Vh);

[P,D] = eig(Iomega_1\Somega);
d = diag(D);
d = sort(d,'descend');
disp('Eigenvalues of S_omega :')
disp(' ');
disp(d(1:10))
dtheo = 1./(2*(0:9)');
dtheo(2) = log(2)/2;
dtheo(1) = 1/2;
err = norm(d(1:10)- dtheo,2)/norm(d(1:10),2);
fprintf('\n Relative error: = %s \n\n',num2str(err));

%% Numerical errors in singular integrals

% Constant rhs : 
T0 = 0*X + 1;
uexact = 2/log(2)*T0;

rhs = integral(Gamma,Vh,T0);
u = Somega\rhs;
figure;
semilogy(m.vtx(:,1),abs(u-uexact(m.vtx)),'--');
title('Inverse of the constant function')
xlabel('x')
ylabel('|u(x) - u_{exact}(x)|')
% Notice that the error near the endpoint of the interval stays quite high.
% Theoretically, the solution being in the Galerkin space, the error should
% be 0 by Céa's lemma. Thus, the error we witness is due to non-exactness
% of the matrix coefficient. To remedy this problem, we would probably need
% an adaptive integration algorithm instead of the choice made here. 

e2 = u'*Iomega_1*u - ...
        2*u'*integral(Gamma,Vh,uexact)...
        + integral(Gamma,uexact^2);
e = sqrt(e2)

%% Order of convergence

Ns = [10; 20; 40; 80; 160; 320; 640];
gss = 3;  

T2 = 2*X^2 - 1; % Chebyshev polynomial
% Now we make it so the solution is not already in the Galerkin space. 
uexact = 4*T2; % Solution of Somega u = T2

e = zeros(length(Ns),1);
for i = 1:length(Ns)
    N = Ns(i);
    m = meshCurve(c,N,'varChange',{@cos,[-pi,0]});
    Gamma = Wdom(m,gss,1/omega,sing);
    Vh = P1(m);
    rhs = integral(Gamma,Vh,T2);
    
    k = 0;
    GXY = @(X,Y)femGreenKernel(X,Y,'[log(r)]',k);
    Somega = -1/(2*pi)*(...
        integral(Gamma,Gamma,Vh,GXY,Vh)  ...
        + regularize(Gamma,Gamma,Vh,'[log(r)]',Vh));
    
    Iomega_1 = integral(Gamma,Vh,Vh);
    
    u = Somega\rhs;
    e2 = u'*Iomega_1*u - ...
        2*u'*integral(Gamma,Vh,uexact)...
        + integral(Gamma,uexact^2);
    e(i) = sqrt(e2);
    
end


figure;

loglog(Ns,e,'-o');
hold on
loglog(Ns,Ns.^(-2),'k--');
legend({'Numerical error','O(1/N^2)'})
xlabel('N')
ylabel('L_{1/\omega}^2 error')
title('Weighted single layer Laplace problem')
eoc = diff(log(e))./diff(log(Ns));
disp('Error:')
disp(e);
disp('Estimated order of convergence')
disp(eoc);

% The error convergence is hindered by the phenomenon described above. 

%% More Gauss points

% If we take more Gauss points, we can verify that this is only due to
% integration error. 

gss = 15;


e = zeros(length(Ns),1);
for i = 1:length(Ns)
    N = Ns(i);
    m = meshCurve(c,N,'varChange',{@cos,[-pi,0]});
    Gamma = Wdom(m,gss,1/omega,sing);
    Vh = P1(m);
    rhs = integral(Gamma,Vh,T2);
    
    k = 0;
    GXY = @(X,Y)femGreenKernel(X,Y,'[log(r)]',k);
    Somega = -1/(2*pi)*(...
        integral(Gamma,Gamma,Vh,GXY,Vh)  ...
        + regularize(Gamma,Gamma,Vh,'[log(r)]',Vh));
    
    Iomega_1 = integral(Gamma,Vh,Vh);
    
    u = Somega\rhs;
    e2 = u'*Iomega_1*u - ...
        2*u'*integral(Gamma,Vh,uexact)...
        + integral(Gamma,uexact^2);
    e(i) = sqrt(e2);
    
end


figure;

loglog(Ns,e,'-o');
hold on
loglog(Ns,Ns.^(-2),'k--');
legend({'Numerical error','O(1/N^2)'})
xlabel('N')
ylabel('L_{1/\omega}^2 error')
title('Weighted single layer Laplace problem')
eoc = diff(log(e))./diff(log(Ns));
disp('Error:')
disp(e);
disp('Estimated order of convergence')
disp(eoc);
