%% Preconditioning of the Symm's integral equation and the hypersingular 
% Laplace equation on the segment [-1,1]



close all;
clear all; %#ok
clc;
tinit = tic;

%% Mesh and boundary element space

c = openline(-1,1);
N = 32000;
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
GMRESTOL = 1e-8;
GMRESMAXIT = 500;
tinit = toc(tinit);

%% 1°) Symm's equation

disp('Symm''s equation');

k = 0;
GXY = @(X,Y)femGreenKernel(X,Y,'[log(r)]',k);

tolEBD = 1e-3;


% Symm's integral operator
if N<5e3
    Somega = -1/(2*pi)*(...
        integral(Gamma,Gamma,Vh,GXY,Vh)  ...
        + regularize(Gamma,Gamma,Vh,'[log(r)]',Vh));
else
    Somega = -1/(2*pi)*(...
        integralEbd(Gamma,Gamma,Vh,'[log(r)]',0,Vh,tolEBD)  ...
        + regularize(Gamma,Gamma,Vh,'[log(r)]',Vh));
end


% Mass matrix int_{Gamma} phi_i phi_j/(1-X^2)^{1/2}
Iomega_1 = integral(Gamma,Vh,Vh);

% Weighted Laplacian
OmegaDx2 = integral(Gamma,grad(Vh),omega2,grad(Vh));

% Galerkin matrix of the L^2_{1/omega} orthogonal projection on constants
L = integral(Gamma,0*X + 1);
C = integral(Gamma,Vh);
pi0 = @(u)((C'*u)/L*C);
% Square-root of the Weighted Laplacian
SQ = @(u)TrefethenSqrtGalerkin(OmegaDx2,5,u,Iomega_1,1,length(Vh)^2);

aux = @(v)(SQ(v) + 2/log(2)*(pi0(v)));
Prec = @(u)(Iomega_1\(aux(Iomega_1\u)));

% rhs "almost singular"
rhs = integral(Gamma,Vh,1/sqrt(1/N^2 + X^2));

tic
[~,~,~,~,RESVEC0]  = gmres(Somega,rhs,[],GMRESTOL,GMRESMAXIT);
toc
tic
[~,~,~,~,RESVEC1]  = gmres(Somega,rhs,[],GMRESTOL,GMRESMAXIT,Iomega_1);
toc
tic
[~,~,~,~,RESVEC2]  = gmres(Somega,rhs,[],GMRESTOL,GMRESMAXIT,Prec);
toc

fprintf('No preconditioner : %s iterations \n',num2str(length(RESVEC0)));
fprintf('Mass matrix preconditioner : %s iterations \n',num2str(length(RESVEC1)));
fprintf('Square-root preconditioner : %s iterations \n',num2str(length(RESVEC2)));

figure; 
semilogy(1:length(RESVEC0),RESVEC0./norm(rhs));
hold on
semilogy(1:length(RESVEC1),RESVEC1./norm(Iomega_1\rhs));
semilogy(1:length(RESVEC2),RESVEC2/norm(Prec(rhs)));
legend({'No preconditioner','Mass matrix','Square-root Laplacian'})
legend boxoff
xlabel('Iteration count')
ylabel('Relative residual')

%% 2°) Hypersingular equation


disp('Hypersingular equation');

% Symm's integral operator
if N < 5e3
    Nomega = -1/(2*pi)*(...
        integral(Gamma,Gamma,omegaDomega(Vh),GXY,omegaDomega(Vh))  ...
        + regularize(Gamma,Gamma,omegaDomega(Vh),'[log(r)]',omegaDomega(Vh)));
else
    Nomega = -1/(2*pi)*(...
        integralEbd(Gamma,Gamma,omegaDomega(Vh),'[log(r)]',0,omegaDomega(Vh),tolEBD)  ...
        + regularize(Gamma,Gamma,omegaDomega(Vh),'[log(r)]',omegaDomega(Vh)));

end

% Mass matrix int_{Gamma} phi_i phi_j/(1-X^2)^{1/2}
Iomega = integral(Gamma,Vh,omega2,Vh);

% Weighted Laplacian
DxOmega2 = integral(Gamma,omegaDomega(Vh),omegaDomega(Vh));

% Galerkin matrix of the L^2_{1/omega} orthogonal projection on constants
L = integral(Gamma,0*X + 1);
C = integral(Gamma,Vh);
% Square-root of the Weighted Laplacian
SQ = @(u)TrefethenSqrtGalerkin(DxOmega2,15,u,Iomega,1.5,length(Vh)^2);


Prec = @(u)(DxOmega2\(SQ(Iomega\u)));

rhs = integral(Gamma,Vh,omega2*(sqrt((1/N^2 + X^2))));

tic
[X0,~,~,~,RESVEC0]  = gmres(Nomega,rhs,[],GMRESTOL,GMRESMAXIT);
toc
tic
[X1,~,~,~,RESVEC1]  = gmres(Nomega,rhs,[],GMRESTOL,GMRESMAXIT,Iomega);
toc
tic
[X2,~,~,~,RESVEC2]  = gmres(Nomega,rhs,[],GMRESTOL,GMRESMAXIT,Prec);
toc

fprintf('No preconditioner : %s iterations \n',num2str(length(RESVEC0)));

fprintf('Mass matrix preconditioner : %s iterations \n',num2str(length(RESVEC1)));
fprintf('Square-root preconditioner : %s iterations \n',num2str(length(RESVEC2)));

figure; 
semilogy(1:length(RESVEC0),RESVEC0./norm(rhs));

hold on
semilogy(1:length(RESVEC1),RESVEC1./norm(Iomega\rhs));
semilogy(1:length(RESVEC2),RESVEC2/norm(Prec(rhs)));

legend({'No preconditioner','Mass matrix','Square-root Laplacian'})
legend boxoff
xlabel('Iteration count')
ylabel('Relative residual')


