% 

% m1 = mshCube(100,[1 1 1]);
m = mshSphere(100,1);
m2 = translate(mshSphere(100,0.2),[1,1,1]);

m1 = union(m,m2);
plot(m1);

Gamma = dom(m,3);
psi_h = fem(m,'P1');
% Wh = fem(m1,'P1');

k = 5;
Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);

M = integral(Gamma,psi_h,psi_h);
dM = integral(Gamma,nx(psi_h),grad(psi_h));

tol = 1e-3;
S = 1/(4*pi)*integral(Gamma,Gamma,psi_h,Gxy,psi_h,tol)...
     + 1/(4*pi)*regularize(Gamma,Gamma,psi_h,'[1/r]',psi_h);

spaceMesh = translate(mshSquare(100,[1,1]),[3,0,0]);

SL = 1/(4*pi)*(integral(spaceMesh.vtx,Gamma,Gxy,psi_h));


hold on
plot(spaceMesh);


PW = @(X)(exp(1i*k*X(:,1)));
rhs = integral(Gamma,psi_h,PW);

lambda = S\rhs;
SL*lambda

plotOn(spaceMesh,real(SL*lambda));