%% Resolution of Helmholtz single layer integral equation on a curve
clear all
close all;
clc;


% c = semicircle; x1 = -1.5; x2 = 1.5; y1 = -1.5; y2 = 1.2; theta_inc = pi/2;
% c = parabola; x1 = -1.5; x2 = 1.5; y1 = -1.5; y2 = 0.4; theta_inc = -pi/6;
c = Scurve; x1 = -1.2; x2 = 1.2; y1 = -1; y2 = 1; theta_inc = -pi/6;
% c = spirale;x1 = -1.5; x2 = 1.5; y1 = -1.5; y2 = 1.8; theta_inc = pi/4;

k = 100;
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

PW = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
rhs = -integral(Gamma,Vh,PW);

GXY = @(X,Y)femGreenKernel(X,Y,'[H0(kr)]',k);
Sh = 1i/4*integral(Gamma,Gamma,Vh,GXY,Vh)  ...
    -1/(2*pi)*regularize(Gamma,Gamma,Vh,'[log(r)]',Vh);


lambda = gmres(Sh,rhs,[],1e-12,size(Sh,1));
figure;
plot(m.vtx(:,1),real(lambda));


x = linspace(x1,x2,2e3);
y = linspace(y1,y2,2e3);
[X,Y] = meshgrid(x,y);
M = [X(:),Y(:),0*X(:)];
NM = size(M,1);

tol = 1e-3;
SL = 1i/4*integralEbd(M,Gamma,'[H0(kr)]',k,Vh,tol)  ...
    -1/(2*pi)*regularize(M,Gamma,'[log(r)]',Vh);

figure,
val = abs(SL*lambda + PW(M));
imagesc(x,y,reshape(val,length(x),length(y)));
axis xy
axis off;
hold on
plot(c);

caxis([0.5,3.5]);

