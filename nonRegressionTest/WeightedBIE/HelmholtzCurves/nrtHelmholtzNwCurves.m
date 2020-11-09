%% Helmholtz single layer integral equation on a segment

clear all; %#ok
close all;
clc




% c = semicircle; x1 = -1.5; x2 = 1.5; y1 = -1.5; y2 = 1.2; theta_inc = pi/2;
% c = parabola; x1 = -1.5; x2 = 1.5; y1 = -1.5; y2 = 0.4; theta_inc = -pi/6;
% c = Scurve; x1 = -1.2; x2 = 1.2; y1 = -1; y2 = 1; theta_inc = -pi/6;
c = spirale;x1 = -1.5; x2 = 1.5; y1 = -1.5; y2 = 1.8; theta_inc = 0;
% c = Vcurve; x1 = -4; x2 = 4; y1 = -2; y2 = 3; theta_inc = pi/2;


k = 50*pi/2; N = fix(5*k*2)+1;
m = meshCurve(c,N,'varChange',{@cos,[-pi,0]});
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
    dOmega{j} = (dw1{j}*w2 + dw2{j}*w1)/(2*omega);
end


singVtx = edges.vtx; % Singularities of 1/omega
singPow = [-1/2;-1/2]; % Power law of the singularities
sing = {singVtx,singPow};

% Integration domain
gss = 3;
Gamma = Wdom(m,gss,1/omega,sing);

Gamma = Gamma.supplyDw(dOmega);

Vh = P1(m);


PW = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
omega2dPW{1} = omega2*1i*k*cos(theta_inc)*PW;
omega2dPW{2} = omega2*1i*k*sin(theta_inc)*PW;
omega2dPW{3} = 0*PW;

rhs = integral(Gamma,ntimes(Vh),omega2dPW);

GXY = @(X,Y)femGreenKernel(X,Y,'[H0(kr)]',k);
omega2GXYomega2 = @(X,Y)(omega2(X).*omega2(Y).*femGreenKernel(X,Y,'[H0(kr)]',k));

N1 = 1i/4*integral(Gamma,Gamma,omegaDomega(Vh),GXY,omegaDomega(Vh))...
    -1/(2*pi)*regularize(Gamma,Gamma,omegaDomega(Vh),'[log(r)]',omegaDomega(Vh));
N2 = -k^2*(...
    1i/4*integral(Gamma,Gamma,ntimes(Vh),omega2GXYomega2,ntimes(Vh))...
    -1/(2*pi)*regularize(Gamma,Gamma,ntimes(Vh),omega2,'omega2[log(r)]',ntimes(Vh)));
Nomega = N1 + N2;
clear N1 N2
% mu = gmres(Nomega,rhs,[],1e-12,size(Nomega,1));
mu = Nomega\rhs;
clear Nomega 
figure;
plot(m.vtx(:,1),real(mu));


x = linspace(x1,x2,1e3);
y = linspace(y1,y2,1e3);
[X,Y] = meshgrid(x,y);
M = [X(:),Y(:),0*X(:)];
NM = size(M,1);

tol = 1e-3;
K = H0Kernel(k);
[Y,Wy] = Gamma.qud;
Ny = size(Y,1);
Wyomega2 = spdiags(Wy.*omega2(Y),0,Ny,Ny);
Mv = uqm(ntimes(Vh),Gamma);
a = 1.5/sqrt(sqrt(NM*Ny));
[Dx,Dy] = offline_dEBD(K,M,Y,a,tol);
A = AbstractMatrix([],...
    @(V)(Dx(Wyomega2*(Mv{1}*V))...
    + Dy(Wyomega2*(Mv{2}*V))),...
    NM,length(Vh));

Reg = regularize(M,Gamma,'grady[log(r)]omega2',ntimes(Vh));
h = figure;
Dmu = 1i/4*(A*mu) - 1/(2*pi)*Reg*mu;
val = abs(Dmu + PW(M));
imagesc(x,y,reshape(val,length(x),length(y)));
axis xy;
axis equal
axis off
hold on
plot(c);
saveas(h,'../spiraleCorlorized');


