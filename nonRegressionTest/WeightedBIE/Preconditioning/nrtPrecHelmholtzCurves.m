%% Square-root preconditioners for the Helmholtz equation on a curved screen

% clearvars;
close all;

%% Mesh and boundary element space

nn = 20;
lambdaEBD = 10;

tolEBD = 1e-3;
GMRESIT = 500;
GMRESTOL = 1e-8;

c = openline(-1,1); theta_inc = pi/4; 
% c = semicircle; x1 = -1.5; x2 = 1.5; y1 = -1.5; y2 = 1.2; theta_inc = pi/2;
% c = parabola; x1 = -1.5; x2 = 1.5; y1 = -1.5; y2 = 0.4; theta_inc = -pi/6;
% c = Scurve; x1 = -1.2; x2 = 1.2; y1 = -1; y2 = 1; theta_inc = -pi/6;
% c = spirale; theta_inc = pi/4;

k = nn/length(c)*pi;
% N = fix(k*length(c)*5)+1;
N = fix(10*k);

GMRESIT = min(GMRESIT,N);


m = meshCurve(c,N,'varChange',{@cos,[-pi,0]});
L = sum(m.ndv);
edges = bnd(m);
% Weight definition :

X1 = edges.vtx(1,:);
X2 = edges.vtx(2,:);
% Weight definition :
[X,Y,Z] = FunR3.XYZ;
w1 = sqrt((X1(1) - X).^2 + (X1(2) - Y).^2 + (X1(3) - Z).^2);
w2 = sqrt((X2(1) - X).^2 + (X2(2) - Y).^2 + (X2(3) - Z).^2);

omega2 = w1*w2;
omega = sqrt(omega2);

dw1{1} = (X - X1(1))./w1;
dw1{2} = (Y - X1(2))./w1;
dw1{3} = (Z - X1(3))./w1;

dw2{1} = (X - X2(1))./w2;
dw2{2} = (Y - X2(2))./w2;
dw2{3} = (Z - X2(3))./w2;

dOmega = cell(1,3);
for j = 1:3
    dOmega{j} = (dw1{j}*w2 + dw2{j}*w1)/(2*sqrt(w1*w2));
end


singVtx = edges.vtx; % Singularities of 1/omega
singPow = [-1/2;-1/2]; % Power law of the singularities
sing = {singVtx,singPow};

% Integration domain
gss = 3;
Gamma = Wdom(m,gss,1/omega,sing);

Gamma = Gamma.supplyDw(dOmega);

Vh = P1(m);


%% 1°) Operator assembling

fprintf('\n\nnn = %s\n',num2str(nn))
fprintf('k = %s\n',num2str(k))
fprintf('N = %s\n',num2str(N))
fprintf('lambda = %s\n',num2str(lambdaEBD))
disp('Assembling operators')

% 
[Somega,Nomega,rq,~] = weightedBIO(Gamma,Vh,k,...
    'lambda',lambdaEBD,'tol',tolEBD,'Ncompress',100);


% GXY = @(X,Y)femGreenKernel(X,Y,'[H0(kr)]',k);
% omega2GXYomega2 = @(X,Y)(omega2(X).*omega2(Y).*GXY(X,Y));
% 
% Nomega = 1i/4*integral(Gamma,Gamma,omegaDomega(Vh),GXY,omegaDomega(Vh))...
%     -1/(2*pi)*regularize(Gamma,Gamma,omegaDomega(Vh),'[log(r)]',omegaDomega(Vh));
% Nomega = Nomega - k^2*(...
%     1i/4*integral(Gamma,Gamma,ntimes(Vh),omega2GXYomega2,ntimes(Vh))...
%     -1/(2*pi)*regularize(Gamma,Gamma,ntimes(Vh),omega2,'omega2[log(r)]',ntimes(Vh)));
% 
% Somega = 1i/4*integral(Gamma,Gamma,Vh,GXY,Vh)  ...
%     -1/(2*pi)*regularize(Gamma,Gamma,Vh,'[log(r)]',Vh);


disp('Done');

% disp('Single layer')

% Somega = 1i/4*integral(Gamma,Gamma,Vh,GXY,Vh)  ...
%     -1/(2*pi)*regularize(Gamma,Gamma,Vh,'[log(r)]',Vh);

omegaDx2 = integral(Gamma,grad(Vh),omega2,grad(Vh));
Omega2 = integral(Gamma,Vh,omega2,Vh);
Iomega_1 = integral(Gamma,Vh,Vh);


K = omegaDx2 - k^2*L^2/4*(Omega2/(L^2/4) - Iomega_1);
keps = k + 1i*0.04*k^(1/3);
SQ1 = @(u)(DarbasPadeSqrt(u,20,pi/3,keps,L^2/4*Iomega_1,K));
PrecSQ1 = @(v)(Iomega_1\SQ1(Iomega_1\v));


dxOmega2 = integral(Gamma,omegaDomega(Vh),omegaDomega(Vh));
Omega2 = integral(Gamma,Vh,omega2^2,Vh);
Iomega = integral(Gamma,Vh,omega2,Vh);

K = dxOmega2 - k^2*L^2/4*(4*Omega2/L^2 - Iomega);
D = dxOmega2 - k^2*Omega2;

SQ2 = @(u)DarbasPadeSqrt(u,15,pi/3,keps,L^2/4*Iomega,K);
PrecSQ2 = @(u)(D\SQ2(Iomega\u));

% Calderon preconditioners

PrecSN = @(u)(Iomega_1\(Somega*(Iomega\u)));
PrecNS = @(u)(Iomega\(Nomega*(Iomega_1\u)));

%% Dirichlet problem:
disp('Dirichlet problem')

PW = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
rhs1 = -integral(Gamma,Vh,PW);


% [Lambda,~,~,~,RESVEC0] = gmres(Somega,rhs1,[],GMRESTOL,GMRESIT);
% [~,~,~,~,RESVEC1] = gmres(Somega,rhs1,[],GMRESTOL,GMRESIT,Iomega_1);
t2 = tic;
[Lambda,~,~,~,RESVEC2] = gmres(Somega,rhs1,[],GMRESTOL,GMRESIT,PrecSQ1);
t2 = toc(t2);
fprintf('Square-root preconditioner : %s iterations in %s seconds \n',...
    num2str(length(RESVEC2)),num2str(t2));

t3 = tic;
[~,~,~,NITER,RESVEC3] = gmres(Somega,rhs1,[],GMRESTOL,GMRESIT,PrecNS);
t3 = toc(t3);
fprintf('Generalized Calderon preconditioner : %s iterations in %s seconds \n',...
    num2str(length(RESVEC3)),num2str(t3));

% 
% figure; 
% % semilogy(1:length(RESVEC0),RESVEC0./norm(rhs1));
% % hold on
% % semilogy(1:length(RESVEC1),RESVEC1./norm(Iomega_1\rhs1));
% semilogy(1:length(RESVEC2),RESVEC2/norm(PrecSQ1(rhs1)));
% semilogy(1:length(RESVEC3),RESVEC3/norm(PrecNS(rhs1)));
% legend({'No preconditioner','Mass matrix',...
%     'Square-root Laplacian','Gen. Calderon Preconditioner'})
% legend boxoff
% xlabel('Iteration count')
% ylabel('Relative residual')

%% 2°) Hypersingular equation


disp('Neumann problem')


PW = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
omega2dPW{1} = omega2*1i*k*cos(theta_inc)*PW;
omega2dPW{2} = omega2*1i*k*sin(theta_inc)*PW;
omega2dPW{3} = 0*PW;

rhs2 = integral(Gamma,ntimes(Vh),omega2dPW);


% [Mu,~,~,~,RESVEC0] = gmres(Nomega,rhs2,[],GMRESTOL,size(Nomega,1));
% [~,~,~,~,RESVEC1] = gmres(Nomega,rhs2,[],GMRESTOL,size(Nomega,1),Iomega);
t2 = tic;
[~,~,~,~,RESVEC2] = gmres(Nomega,rhs2,[],GMRESTOL,GMRESIT,PrecSQ2);
t2 = toc(t2);
fprintf('Square-root preconditioner : %s iterations in %s seconds \n',...
    num2str(length(RESVEC2)),num2str(t2));

t3 = tic;
[~,~,~,~,RESVEC3] = gmres(Nomega,rhs2,[],GMRESTOL,GMRESIT,PrecSN);
t3 = toc(t3);
fprintf('Generalized Calderon preconditioner : %s iterations in %s seconds \n',...
    num2str(length(RESVEC3)),num2str(t3));


% fprintf('No preconditioner : %s iterations \n',num2str(length(RESVEC0)));
% fprintf('Mass matrix preconditioner : %s iterations \n',num2str(length(RESVEC1)));
% fprintf('Square-root preconditioner : %s iterations \n',num2str(length(RESVEC2)));
% fprintf('Generalized Calderon preconditioner : %s iterations \n',num2str(length(RESVEC3)));

% figure; 
% % semilogy(1:length(RESVEC0),RESVEC0/norm(rhs2));
% hold on
% % semilogy(1:length(RESVEC1),RESVEC1/norm((Iomega\rhs2)));
% semilogy(1:length(RESVEC2),RESVEC2/norm(PrecSQ2(rhs2)));
% semilogy(1:length(RESVEC3),RESVEC3/norm(PrecSN(rhs2)));
% legend({'No preconditioner','Mass matrix','Square-root Laplacian','Gen. Calderon Preconditioner'})
% legend boxoff
% xlabel('Iteration count')
% ylabel('Relative residual')

