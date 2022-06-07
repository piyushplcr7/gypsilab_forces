addpath(genpath("../../"));
clear; 
clc;
format long;

N = 200;
rad = 2;
mesh = mshSphere(N,rad);
order = 3;
Gamma = dom(mesh,order);

% parameters
stp = mesh.stp;
% Wave number, stable?
k = 1/stp(2);
f = (k*340)/(2*pi);
omega = k;

% FEM spaces
S0_Gamma = fem(mesh,'P0');
S1_Gamma = fem(mesh,'P1');
E0x_Gamma = fem(mesh,'RWG');

% Kernel
Gxy = @(X,Y)femGreenKernel(X,Y,'[exp(ikr)/r]',k);

% Constructing the block matrix

% Top left
A = integral(Gamma,Gamma,E0x_Gamma,Gxy,E0x_Gamma);
A = A + regularize(Gamma,Gamma,E0x_Gamma,'[1/r]',E0x_Gamma);
A = A * (-1./4/pi);

% Top right
B = integral(Gamma,Gamma,div(E0x_Gamma),Gxy,S0_Gamma);
B = B + regularize(Gamma,Gamma,div(E0x_Gamma),'[1/r]',S0_Gamma);
B = B * (-1./4/pi);

% Bottom left
C = B';

% Bottom right identity
D = integral(Gamma,Gamma,S0_Gamma,Gxy,S0_Gamma);
D = D + regularize(Gamma,Gamma,S0_Gamma,'[1/r]',S0_Gamma);
D = D * (-1./4/pi) * k^2;

% Block matrix
T = [A B;
     C D];

% Constructing the RHS

% Plane waves
% Incident direction and field
X0 = [0 0 -1]; 
E0  = [0 1  0]; % Polarization (+x for Theta-Theta and +y for Phi-Phi)
% H0  = cross(X0,E);

% % Incident Plane wave (electromagnetic field) which is used as the
% % potential U
% Uin{1} = @(X) 1i/omega * exp(1i*k*X*X0') * E0(1);
% Uin{2} = @(X) 1i/omega * exp(1i*k*X*X0') * E0(2);
% Uin{3} = @(X) 1i/omega * exp(1i*k*X*X0') * E0(3);

% Incident Plane wave (electromagnetic field) which is used as the
% potential U
Ein{1} = @(X) exp(1i*k*X*X0') * E0(1);
Ein{2} = @(X) exp(1i*k*X*X0') * E0(2);
Ein{3} = @(X) exp(1i*k*X*X0') * E0(3);

% Incident currents
nxEin = - integral(Gamma,E0x_Gamma,E0x_Gamma)\(integral(Gamma,nx(E0x_Gamma),Ein));
nxUin = 1i/omega * nxEin;

% Making the assumption (needs to be verified) that our
% L1(n)=(((1/2)Id+C) (nx(Uin)xn),v) is equal to (((1/2)Id + M)nxUin,nxv) !!!!!

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
gradxGxy{1} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]1',k) ;
gradxGxy{2} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]2',k) ;
gradxGxy{3} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]3',k) ;

M = - 1/(4*pi) * integral(Gamma, Gamma,nx(E0x_Gamma), gradxGxy, E0x_Gamma); 
M = M + 1/(4*pi) * regularize(Gamma, Gamma,nx(E0x_Gamma),'grady[1/r]', E0x_Gamma);

Id = integral(Gamma,E0x_Gamma,E0x_Gamma);

% Solving integral equation
L1g = (1/2*Id + M)*nxUin;

rhsvec = [L1g; zeros(size(S0_Gamma.dof,1),1)];

sol = T\rhsvec;

%udotn = sol(size(E0x_Gamma.dof,1)+1:size(sol));
j = -sol(1:size(E0x_Gamma.dof,1));


% RESULTING TOTAL ELECTRIC FIELD
% Use Stratton-Chu representation

radiat = mshSquare(1e3,[4 4]);

% Esc = mathcal{I}j - matchal{K}m
% where j=-nxHsc is the solution above and m = nxEinc
v = E0x_Gamma;


I1 = 1i * k * 1/(4*pi)*integral(radiat.vtx, Gamma, Gxy,rest(v,1))... 
        + 1i/k * 1/(4*pi) * integral(radiat.vtx, Gamma,gradxGxy{1}, div(v));

I2 = 1i * k * 1/(4*pi)*integral(radiat.vtx, Gamma, Gxy,rest(v,2))... 
        + 1i/k * 1/(4*pi) * integral(radiat.vtx, Gamma, gradxGxy{2}, div(v));

I3 = 1i * k * 1/(4*pi)*integral(radiat.vtx, Gamma, Gxy,rest(v,3))... 
        + 1i/k * 1/(4*pi) * integral(radiat.vtx, Gamma, gradxGxy{3}, div(v));

    
K1 = 1/(4*pi) * integral(radiat.vtx, Gamma, gradxGxy{2},rest(v,3))...
        - 1/(4*pi) * integral(radiat.vtx, Gamma, gradxGxy{3},rest(v,2));
K2 = 1/(4*pi) * integral(radiat.vtx, Gamma, gradxGxy{3},rest(v,1))...
        - 1/(4*pi) * integral(radiat.vtx, Gamma, gradxGxy{1},rest(v,3));
K3 = 1/(4*pi) * integral(radiat.vtx, Gamma, gradxGxy{1},rest(v,2))...
        - 1/(4*pi) * integral(radiat.vtx, Gamma, gradxGxy{2},rest(v,1));
    
Esc1 = I1*j - K1*nxEin;
Esc2 = I2*j - K2*nxEin;
Esc3 = I3*j - K3*nxEin;

EtestNorm = sqrt(real(Esc1).^2 + real(Esc2).^2 + real(Esc3).^2);
figure

vtx = radiat.vtx;
X = vtx(:,1); Y = vtx(:,2); Z = vtx(:,3);
U = real(Esc1); V = real(Esc2);  W = real(Esc3); 
quiver3(X,Y,Z,U,V,W,'r')
hold on
radiat.vtx(:,3) = EtestNorm;
plot(radiat,EtestNorm)
axis equal
title('|j|')
xlabel('X');   ylabel('Y');   zlabel('Z');
colorbar


% test space first, trial space second
%size(integral(Gamma,Gamma,S0_Gamma,Gxy,S1_Gamma))

