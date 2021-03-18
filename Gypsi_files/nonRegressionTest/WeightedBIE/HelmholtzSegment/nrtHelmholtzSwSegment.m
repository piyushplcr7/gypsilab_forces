%% Helmholtz single layer integral equation on a segment

clear all
close all;
clc

k = 100;
c = unitSegment;
m = meshCurve(c,500,'varChange',{@cos,[-pi,0]});
edges = bnd(m);
X1 = edges.vtx(1,:);
X2 = edges.vtx(2,:);
% Weight definition :
[X,Y,Z] = FunR3.XYZ;
omega1 = sqrt((X1(1) - X).^2 + (X1(2) - Y).^2 + (X1(3) - Z).^2);
omega2 = sqrt((X2(1) - X).^2 + (X2(2) - Y).^2 + (X2(3) - Z).^2);
omega = sqrt(omega1*omega2);
singVtx = [X1;X2]; % Singularities of 1/omega
singPow = [-1/2;-1/2]; % Power law of the singularities
sing = {singVtx,singPow};

% Integration domain
gss = 5;
Gamma = Wdom(m,gss,1/omega,sing);

Vh = P1(m);

theta_inc = -pi/6;

PW = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
rhs = -integral(Gamma,Vh,PW);

GXY = @(X,Y)femGreenKernel(X,Y,'[H0(kr)]',k);
Sh = 1i/4*integral(Gamma,Gamma,Vh,GXY,Vh)  ...
    -1/(2*pi)*regularize(Gamma,Gamma,Vh,'[log(r)]',Vh);


lambda = gmres(Sh,rhs,[],1e-12,size(Sh,1));
figure;
plot(m.vtx(:,1),real(lambda));

radiat = mshSquare(5e5,[6,4]);

tol = 1e-3;
SLnonreg = integralEbd(radiat.vtx,Gamma,'[H0(kr)]',k,Vh,tol);
SLreg = regularize(radiat.vtx,Gamma,'[log(r)]',Vh);

SL = 1i/4*SLnonreg + (-1)/(2*pi)*SLreg;

figure,
m2 = radiat;
m2.vtx(:,3) = abs(SL*lambda + PW(radiat.vtx));
plotOn(m2,m2.vtx(:,3));
axis equal
axis off
caxis([0.5,3]);

