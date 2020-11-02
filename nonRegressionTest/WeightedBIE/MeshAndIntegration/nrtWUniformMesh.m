% Creation of uniform meshes with respect to a weight function
% Computation of Gaussian quadrature adapted to the weight function. 


%% Chebyshev weight 

m = openline(-1,1);
m = meshCurve(m,100);
[~,B] = m.ABCD; % 
singVtx = [-1,0,0;1,0,0];
singPow = [-1/2;-1/2];
[X,Y,Z] = FunR3.XYZ;
omega = 1/sqrt(1-X^2);
p = WprimitiveOnMesh(m,B,omega,singVtx,singPow,5);
sum(p) - pi
plot(pi - acos(m.vtx(:,1)));
hold on
plot(cumsum([0;p]),'--');


%% Uniform mesh adapted to the weight

clear all 
close all

c = openline(-1,1);
singVtx = [1 0 0; -1 0 0; 0 0 0];
singPow = [-1/2;-1/2;-0.8];
[X,Y,Z] = FunR3.XYZ;
omega = 1/sqrt(1 - X^2).*(abs(X)).^(singPow(3));
m = WuniformMesh(c,50,omega,singVtx,singPow);

plot(m);


%% Weight-uniform mesh of a spirale-shaped arc.  

c = spirale;
c = c.normalParam;
m = meshCurve(c,100);
boundary = bnd(m);
[~,B] = m.ABCD;
singVtx = boundary.vtx;
X1 = singVtx(1,:);
X2 = singVtx(2,:);
dX1 = sqrt((X-X1(1))^2 + (Y - X1(2))^2);
dX2 = sqrt((X-X2(1))^2 + (Y - X2(2))^2);
omega = 1/sqrt(dX1*dX2);
singPow = [-1/2;-1/2];

p = WprimitiveOnMesh(m,B,omega,singVtx,singPow,20);
plot(cumsum([0;p]));

m = WuniformMesh(c,200,omega,singVtx,singPow);
[~,B] = m.ABCD;
p = WprimitiveOnMesh(m,B,omega,singVtx,singPow,10);

figure
plot(m);
axis equal

%% Computation of Gaussian quadrature for the weight


c = unitSegment;
m = meshCurve(c,100);
boundary = bnd(m);
[~,B] = m.ABCD;
singVtx = [boundary.vtx];
singPow = [-1/2;-1/2];
omega = 1/(sqrt(1-X^2));

m = WuniformMesh(c,20,omega,singVtx,singPow);
[~,B] = m.ABCD;
figure; 
plot(m);
axis equal

gss = 3;
sing = {singVtx,singPow};
Gamma = Wdom(m,gss,omega,sing);

[X,Y,Z] = FunR3.XYZ;
% int_{-1}^1 sin(x)/sqrt(1 - x^2)
int_approx = integral(Gamma,sin(X).^2);
% Reference value compute by using x = cos(t) change of variables.
int_ref = integral(@(t)(sin(cos(t)).^2),-pi,0);
int_approx - int_ref

