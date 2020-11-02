%% Resolution of Symm's integral equation
clear all
close all;
clc;


%% Mesh and integration domain

c = openline(-1,1);
m = meshCurve(c,400,'varChange',{@cos,[-pi,0]});
edges = bnd(m);
% Weight definition :
[X,Y,Z] = FunR3.XYZ;
omega2 = 1 - X^2;
omega = sqrt(omega2);

dOmega{1} = -X./omega;
dOmega{2} = 0*X;
dOmega{3} = 0*X;

singVtx = edges.vtx; % Singularities of 1/omega
singPow = [-1/2;-1/2]; % Power law of the singularities
sing = {singVtx,singPow};

% Integration domain
gss = 5;
Gamma = Wdom(m,gss,1/omega,sing);
Gamma = Gamma.supplyDw(dOmega);
check = integral(Gamma,omega2);
pi/2 - check

% Boundary element space
Vh = P1(m);

%% Generalized eigenvalue problem


k = 0;
GXY = @(X,Y)femGreenKernel(X,Y,'[log(r)]',k);

omega2GXYomega2 = @(X,Y)(omega2(X).*omega2(Y).*femGreenKernel(X,Y,'[log(r)]',k));
GXYomega2 = @(X,Y)(omega2(Y).*femGreenKernel(X,Y,'[log(r)]',k));
omega2GXY = @(X,Y)(omega2(X).*femGreenKernel(X,Y,'[log(r)]',k));


% Alternative way : 
Nomega = -1/(2*pi)*(integral(Gamma,Gamma,omegaDomega(Vh),GXY,omegaDomega(Vh))...
    + regularize(Gamma,Gamma,omegaDomega(Vh),'[log(r)]',omegaDomega(Vh)));

Iomega = integral(Gamma,Vh,omega2,Vh);


[P,D] = eig(Iomega\Nomega);
d = diag(D);
d = sort(d,'ascend');
disp(d(1:10));
dtheo = (1:10)'/2;
err = norm(d(1:10)- dtheo,2)/norm(d(1:10),2);

fprintf('\n Relative error: = %s \n\n',num2str(err));

%% Numerical errors in singular integrals

U0 = 0*X + 1;
uexact = 2*U0; % Solution of Nomega u = U0

rhs = integral(Gamma,Vh,omega2*U0);
u = Nomega\rhs;

figure;
semilogy(m.vtx(:,1),abs(u-uexact(m.vtx)),'--');
title('Inverse of the constant function')
xlabel('x')
ylabel('|u(x) - u_{exact}(x)|')
% Notice that the error near the endpoint of the interval stays quite high.
% Theoretically, the solution being in the Galerkin space, the error should
% be 0 by Céa's lemma. Thus, the error we witness is due to non-exactness
% of the matrix coefficient.
% It is not a big problem in the Hypersingular setup since the weight omega
% is small near the endpoints.

e2 = u'*Iomega*u - ...
    2*u'*integral(Gamma,Vh,omega2*uexact)...
    + integral(Gamma,omega2*uexact^2);

e = sqrt(e2)

%% Order of convergence

Ns = [10; 20; 40; 80; 160; 320; 640];

gss = 3;

U2 = 4*X^2 - 1; % Chebyshev polynomial of second kind
% Now we make it so the solution is not already in the Galerkin space.
uexact = 2/3*U2; % Solution of Nomega u = U2

e = zeros(length(Ns),1);
for i = 1:length(Ns)
    
    N = Ns(i);
    m = meshCurve(c,N,'varChange',{@cos,[-pi,0]});
    Gamma = Wdom(m,gss,1/omega,sing);
    Gamma = Gamma.supplyDw(dOmega);
    Vh = P1(m);
    rhs = integral(Gamma,Vh,omega2*U2);
    
    k = 0;
    Nomega = -1/(2*pi)*(integral(Gamma,Gamma,omegaDomega(Vh),GXY,omegaDomega(Vh))...
    + regularize(Gamma,Gamma,omegaDomega(Vh),'[log(r)]',omegaDomega(Vh)));
    
    Iomega = integral(Gamma,Vh,omega2,Vh);
    
    
    u = Nomega\rhs;
    e2 = u'*Iomega*u - ...
        2*u'*integral(Gamma,Vh,omega2*uexact)...
        + integral(Gamma,omega2*uexact^2);
    e(i) = sqrt(e2);
    
end


figure;

loglog(Ns,e,'-o');
hold on
loglog(Ns,Ns.^(-2),'k--');
legend({'Numerical error','O(1/N^2)'})
xlabel('N')
ylabel('L_{\omega}^2 error')
title('Weighted hypersingular Laplace problem')
eoc = diff(log(e))./diff(log(Ns));
disp('Error:')
disp(e);
disp('Estimated order of convergence')
disp(eoc);


