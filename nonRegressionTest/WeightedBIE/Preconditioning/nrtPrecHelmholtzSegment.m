%% Square-root preconditioners for the Helmholtz equation on a flat screen. 

clear all;
close all;
clc;

GMRESTOL = 1e-8;

%% Mesh and boundary element space

k = 200*pi/2; N = fix(5*k*2);

c = openline(-1,1);
m = meshCurve(c,N,'varChange',{@cos,[-pi,0]});
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
gss = 3;
Gamma = Wdom(m,gss,1/omega,sing);
Gamma = Gamma.supplyDw(dOmega);


% Boundary element space
Vh = P1(m);

%% 1°) Single-layer equation

disp('Single layer integral equation')

theta_inc = 0;

PW = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
rhs1 = -integral(Gamma,Vh,PW);

GXY = @(X,Y)femGreenKernel(X,Y,'[H0(kr)]',k);
Somega = 1i/4*integral(Gamma,Gamma,Vh,GXY,Vh)  ...
    -1/(2*pi)*regularize(Gamma,Gamma,Vh,'[log(r)]',Vh);


omegaDx2 = integral(Gamma,grad(Vh),omega2,grad(Vh));
Omega2 = integral(Gamma,Vh,omega2,Vh);
Iomega_1 = integral(Gamma,Vh,Vh);

K = omegaDx2 - k^2*(Omega2 - Iomega_1);
keps = k + 1i*0.025*k^(1/3);
SQ = @(u)(DarbasPadeSqrt(u,20,pi/3,keps,Iomega_1,K));
Prec = @(v)(Iomega_1\SQ(Iomega_1\v));

[~,~,~,~,RESVEC0] = gmres(Somega,rhs1,[],GMRESTOL,size(Somega,1));
[~,~,~,~,RESVEC1] = gmres(Iomega_1\Somega,Iomega_1\rhs1,[],GMRESTOL,size(Somega,1));
[~,~,~,~,RESVEC2] = gmres(Somega,rhs1,[],GMRESTOL,size(Somega,1),Prec);


fprintf('No preconditioner : %s iterations \n',num2str(length(RESVEC0)));
fprintf('Mass matrix preconditioner : %s iterations \n',num2str(length(RESVEC1)));
fprintf('Square-root preconditioner : %s iterations \n',num2str(length(RESVEC2)));

figure; 
semilogy(1:length(RESVEC0),RESVEC0./norm(rhs1));
hold on
semilogy(1:length(RESVEC1),RESVEC1./norm(Iomega_1\rhs1));
semilogy(1:length(RESVEC2),RESVEC2/norm(Prec(rhs1)));
legend({'No preconditioner','Mass matrix','Square-root Laplacian'})
legend boxoff
xlabel('Iteration count')
ylabel('Relative residual')

%% 2°) Hypersingular equation


disp('Hypersingular integral equation');

theta_inc = -pi/4;

PW = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
omega2dPW{1} = omega2*1i*k*cos(theta_inc)*PW;
omega2dPW{2} = omega2*1i*k*sin(theta_inc)*PW;
omega2dPW{3} = 0*PW;

rhs2 = integral(Gamma,ntimes(Vh),omega2dPW);

GXY = @(X,Y)femGreenKernel(X,Y,'[H0(kr)]',k);
omega2GXYomega2 = @(X,Y)(omega2(X).*omega2(Y).*femGreenKernel(X,Y,'[H0(kr)]',k));

N1 = 1i/4*integral(Gamma,Gamma,omegaDomega(Vh),GXY,omegaDomega(Vh))...
    -1/(2*pi)*regularize(Gamma,Gamma,omegaDomega(Vh),'[log(r)]',omegaDomega(Vh));
N2 = -k^2*(...
    1i/4*integral(Gamma,Gamma,ntimes(Vh),omega2GXYomega2,ntimes(Vh))...
    -1/(2*pi)*regularize(Gamma,Gamma,ntimes(Vh),omega2,'omega2[log(r)]',ntimes(Vh)));
Nomega = N1 + N2;

dxOmega2 = integral(Gamma,omegaDomega(Vh),omegaDomega(Vh));
Omega2 = integral(Gamma,Vh,omega2^2,Vh);
Iomega = integral(Gamma,Vh,omega2,Vh);

K = dxOmega2 - k^2*(Omega2 - Iomega);
D = dxOmega2 - k^2*Omega2;

SQ = @(u)DarbasPadeSqrt(u,15,pi/3,keps,Iomega,K);

Prec = @(u)(D\SQ(Iomega\u));

[~,~,~,~,RESVEC0] = gmres(Nomega,rhs2,[],1e-12,size(Nomega,1));
[~,~,~,~,RESVEC1] = gmres(Nomega,rhs2,[],1e-12,size(Nomega,1),Iomega);
[~,~,~,~,RESVEC2] = gmres(Nomega,rhs2,[],1e-12,size(Nomega,1),Prec);



fprintf('No preconditioner : %s iterations \n',num2str(length(RESVEC0)));
fprintf('Mass matrix preconditioner : %s iterations \n',num2str(length(RESVEC1)));
fprintf('Square-root preconditioner : %s iterations \n',num2str(length(RESVEC2)));

figure; 
semilogy(1:length(RESVEC0),RESVEC0);
hold on
semilogy(1:length(RESVEC1),RESVEC1/norm((Iomega\rhs2)));
semilogy(1:length(RESVEC2),RESVEC2/norm(Prec(rhs2)));
legend({'No preconditioner','Mass matrix','Square-root Laplacian'})
legend boxoff
xlabel('Iteration count')
ylabel('Relative residual')

%% 3°) Calderon preconditioning

disp('Calderon preconditioning')

PrecS = @(u)(Iomega_1\(Somega*(Iomega\u)));
PrecN = @(u)(Iomega\(Nomega*(Iomega_1\u)));


[~,~,~,~,RESVEC0] = gmres(Somega,rhs1,[],1e-12,size(Somega,1),PrecN);
[~,~,~,~,RESVEC1] = gmres(Nomega,rhs2,[],1e-12,size(Nomega,1),PrecS);

fprintf('NS preconditioning: %s iterations \n',num2str(length(RESVEC0)));
fprintf('SN preconditioning: %s iterations \n',num2str(length(RESVEC1)));



