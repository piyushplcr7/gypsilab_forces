%
clear all
close all
clc;

[X,Y,Z] = FunR3.XYZ;
m = bnd(mshCube(800,[1 1 1]));
m = rotate(m,[1 1 1],pi/4);
% m = mshTorus(800,1,0.5);
wall = rotate(mshSquare(800,[4,4]),[1 0 0],pi/2);
wall = translate(wall,[0 2 0]);

k = 12;
Gk = @(X,Y)(femGreenKernel(X,Y,'[exp(ikr)/r]',k));
G0 = @(X,Y)(femGreenKernel(X,Y,'[1/r]',0));
dGk = cell(3,1);
dGk{1} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k));
dGk{2} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k));
dGk{3} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k));

Vh = fem(m,'P1');
Gamma = dom(m,3);

Sk = 1/(4*pi)*integral(Gamma,Gamma,Vh,Gk,Vh) ...
    + 1/(4*pi)*regularize(Gamma,Gamma,Vh,'[1/r]',Vh);

Dk = 1/(4*pi)*integral(Gamma,Gamma,Vh,dGk,ntimes(Vh))...
    + 1/(4*pi)*regularize(Gamma,Gamma,Vh,'grady[1/r]',ntimes(Vh));


Nk = 1/(4*pi)*integral(Gamma,Gamma,nxgrad(Vh),Gk,nxgrad(Vh)) ...
    + 1/(4*pi)*regularize(Gamma,Gamma,nxgrad(Vh),'[1/r]',nxgrad(Vh))...
    -k^2/(4*pi)*integral(Gamma,Gamma,ntimes(Vh),Gk,ntimes(Vh))...
    -k^2/(4*pi)*regularize(Gamma,Gamma,ntimes(Vh),'[1/r]',ntimes(Vh));

N0 = 1/(4*pi)*integral(Gamma,Gamma,nxgrad(Vh),G0,nxgrad(Vh)) ...
    + 1/(4*pi)*regularize(Gamma,Gamma,nxgrad(Vh),'[1/r]',nxgrad(Vh));



M = integral(Gamma,Vh,Vh);
Delta = integral(Gamma,grad(Vh),grad(Vh));

keps = k + 1i*0.05*k^(1/3);
DtNk = -padePrecondDarbas([],20,pi/2,keps,M,Delta);
I = eye(size(M,1));
CFIE = (M\Sk)*(M\DtNk) + (I/2 + M\Dk);


NoPrec = eye(size(Sk,1));
Prec0 = M^(-1)*(N0+M)*M^(-1);
Preck = M^(-1)*(Nk)*M^(-1);

PW = exp(1i*k*(Y+Z)/sqrt(2));

rhs = integral(Gamma,Vh,PW);

[u,~,~,~,RESVEC1] = gmres(NoPrec*Sk,NoPrec*rhs,[],1e-8,500);
[~,~,~,~,RESVEC2] = gmres(Prec0*Sk,Prec0*rhs,[],1e-8,500);
[~,~,~,~,RESVEC3] = gmres(Preck*Sk,Preck*rhs,[],1e-8,500);
[mu,~,~,~,RESVEC4] = gmres(CFIE,M\rhs,[],1e-8,500);


% Sanity check CFIE:
SL = 1/(4*pi)*integral(wall.vtx,Gamma,Gk,Vh) ...
    + 1/(4*pi)*regularize(wall.vtx,Gamma,'[1/r]',Vh);

DL = 1/(4*pi)*integral(wall.vtx,Gamma,dGk,ntimes(Vh)) ...
    + 1/(4*pi)*regularize(wall.vtx,Gamma,'grady[1/r]',ntimes(Vh));

figure;
plotOn(wall,real(SL*u - PW(wall.vtx)));
view(0,0);
hold on
plotOn(m,real((M\Sk)*u - PW(m.vtx)));
p1 = [0 0 0];                         % First Point
p2 = [0 1 1];                         % Second Point
dp = p2-p1;                         % Difference
quiver3(p1(1),p1(2),dp(1),dp(2),0);
figure
plotOn(wall,real(SL*(M\DtNk)*mu + DL*mu - PW(wall.vtx)));
view(0,0);
hold on
plotOn(m,real((M\Sk)*u - PW(m.vtx)));


figure
semilogy(RESVEC1/norm(NoPrec*rhs));
hold on
semilogy(RESVEC2/norm(Prec0*rhs));
semilogy(RESVEC3/norm(Preck*rhs));
semilogy(RESVEC4/norm(M\rhs));

figure
plotOn(m,real(u));