%% Resolution of Helmholtz single layer integral equation on a curve
clear all
close all;
clc;


% c = semicircle; x1 = -1.5; x2 = 1.5; y1 = -1.5; y2 = 1.2; theta_inc = pi/2;
% c = parabola; x1 = -1.5; x2 = 1.5; y1 = -1.5; y2 = 0.4; theta_inc = -pi/6;
% c = Scurve; x1 = -1.5; x2 = 1.5; y1 = -0.5; y2 = 2; theta_inc = -pi/2;
% c = spirale;x1 = -1.5; x2 = 1.5; y1 = -1.5; y2 = 1.8; theta_inc = pi/4;
c = Vcurve; x1 = -2; x2 = 2; y1 = -1; y2 = 3; theta_inc = pi/2;
% L = length(c);
L=2;
% X1 = [c.x(c.I(1)),c.y(c.I(1)),0];
% X2 = [c.x(c.I(2)),c.y(c.I(2)),0];
% X3 = [c.x(0),c.y(0),0];
% [X,Y,Z] = FunR3.XYZ;
% omega1 = sqrt((X1(1) - X).^2 + (X1(2) - Y).^2 + (X1(3) - Z).^2);
% omega2 = sqrt((X2(1) - X).^2 + (X2(2) - Y).^2 + (X2(3) - Z).^2);
% omega3 = sqrt((X3(1) - X).^2 + (X3(2) - Y).^2 + (X3(3) - Z).^2);
% omega = sqrt(omega1*omega2)*omega3^(1/3);
% singVtx = [X1;X2;X3];
% singPow = [-1/2;-1/2; -1/3];

c = Vcurve; x1 = -2; x2 = 2; y1 = -1; y2 = 3; theta_inc = pi/2;
L = length(c);
L=2;
X1 = [c.x(c.I(1)),c.y(c.I(1)),0];
X2 = [c.x(c.I(2)),c.y(c.I(2)),0];
[X,Y,Z] = FunR3.XYZ;
omega1 = sqrt((X1(1) - X).^2 + (X1(2) - Y).^2 + (X1(3) - Z).^2);
omega2 = sqrt((X2(1) - X).^2 + (X2(2) - Y).^2 + (X2(3) - Z).^2);
omega = sqrt(omega1*omega2);
singVtx = [X1;X2];
singPow = [-1/2;-1/2];

k = 50*pi/L;
N = fix(5*k*L)+2;
m = WuniformMesh(c,N,1/omega,singVtx,singPow);

sing{1} = singVtx;
sing{2} = singPow;

% Integration domain
gss = 3;
Gamma = Wdom(m,gss,1/omega,sing);

Vh = P1(m);

PW = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
rhs = -integral(Gamma,Vh,PW);

GXY = @(X,Y)femGreenKernel(X,Y,'[H0(kr)]',k);
Sh = 1i/4*integral(Gamma,Gamma,Vh,GXY,Vh)  ...
    -1/(2*pi)*regularize(Gamma,Gamma,Vh,'[log(r)]',Vh);


omegaDx2 = integral(Gamma,grad(Vh),omega^2,grad(Vh));
Omega2 = integral(Gamma,Vh,omega^2,Vh);
Iomega_1 = integral(Gamma,Vh,Vh);


K = omegaDx2 - k^2*L^2/4*(Omega2/(L^2/4) - Iomega_1);
keps = k + 1i*0.04*k^(1/3);
SQ1 = @(u)(DarbasPadeSqrt(u,20,pi/3,keps,L^2/4*Iomega_1,K));
PrecSQ1 = @(v)(Iomega_1\SQ1(Iomega_1\v));

lambda = gmres(Sh,rhs,[],1e-8,500,PrecSQ1);
plot(m.vtx(:,1),real(lambda));


x = linspace(x1,x2,2e3);
y = linspace(y1,y2,2e3);
[X,Y] = meshgrid(x,y);
M = [X(:),Y(:),0*X(:)];
NM = size(M,1);

tol = 1e-3;
lambdaEBD = 5;
SL = 1i/4*integralEbd(M,Gamma,'[H0(kr)]',k,Vh,tol,lambdaEBD)  ...
    -1/(2*pi)*regularize(M,Gamma,'[log(r)]',Vh);

figure;
val = abs(SL*(lambda) + PW(M));
imagesc(x,y,reshape(val,length(x),length(y)));
axis xy
axis off;
hold on
plot(c);

caxis([.2,2.5]);

