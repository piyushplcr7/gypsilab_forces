%+========================================================================+
% Maxwell electric wave equation by EFIE
% (Extension by -Einc approach in the interior)
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
%run('../addpathGypsilab.m')
addpath(genpath("../../"));

% Parameters
N   = 1e3; 
tol = 1e-5; 
gss = 3; 

% Spherical mesh
mesh = mshSphere(N,1);
Gamma  = dom(mesh,3);
% figure
% plot(mesh)
% hold on
% plotNrm(mesh)

% Finite elements
v = fem(mesh,'RWG');

% Setting a frequency suitable for testing
stp = mesh.stp;
k   = 3;
c   = 299792458;
f   = (k*c)/(2*pi);
disp(['Frequency : ',num2str(f/1e6),' MHz']);

% Incident direction and field
X0 = [0 0 -1]; 
E  = [0 1  0]; 
H  = cross(X0,E);

% Incident plane wave (polarized electromagnetic field)
PWE{1} = @(X) exp(1i*k*X*X0') * E(1);
PWE{2} = @(X) exp(1i*k*X*X0') * E(2);
PWE{3} = @(X) exp(1i*k*X*X0') * E(3);

PWH{1} = @(X) exp(1i*k*X*X0') * H(1);
PWH{2} = @(X) exp(1i*k*X*X0') * H(2);
PWH{3} = @(X) exp(1i*k*X*X0') * H(3);

% Plot incident wave
figure
plot(mesh)
axis equal
plot(mesh,real(PWE{2}(mesh.vtx)))
title('Incident wave')
xlabel('X');   ylabel('Y');   zlabel('Z');
view(0,10)
hold off


%%% PREPARE OPERATORS
disp('~~~~~~~~~~~~~ PREPARE OPERATOR ~~~~~~~~~~~~~')

disp('Computing convolution kernels')
% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
Gxy = @(X,Y) femGreenKernel(X,Y,'[exp(ikr)/r]',k);

% Green kernel function --> G(x,y) = exp(ik|x-y|)/|x-y| 
gradxGxy{1} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]1',k) ;
gradxGxy{2} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]2',k) ;
gradxGxy{3} = @(X,Y) femGreenKernel(X,Y,'gradx[exp(ikr)/r]3',k) ;


% Finite element boundary operator
disp('Assembling T')
tic
T = k^2/(4*pi)*integral(Gamma, Gamma, v, Gxy, v, tol) ...
    - 1/(4*pi)*integral(Gamma, Gamma, div(v), Gxy, div(v), tol) ;
toc

disp('Regularizing T')
tic
T = T + k^2/(4*pi)*regularize(Gamma, Gamma, v, '[1/r]', v) ...
      - 1/(4*pi)*regularize(Gamma, Gamma, div(v), '[1/r]', div(v));
toc

Id = integral(Gamma,v,v);

% Incident currents 
nxEinc = - integral(Gamma,v,v)\(integral(Gamma,nx(v),PWE)); 
nxHinc = - integral(Gamma,v,v)\(integral(Gamma,nx(v),PWH)); 

RHS = -integral(Gamma,v,PWE);

disp('Solving integral equation')
j = -1i*k*(T \ RHS);


%%% SURFACE CURRENT
disp('~~~~~~~~~~~~~ Surface current ~~~~~~~~~~~~~')

% Mesh Interpolation
jmsh = feval(v, j,mesh);
jnorm    = sqrt(sum(real(cell2mat(jmsh)).^2,2));

% Graphical representation
figure

vtx = mesh.vtx;
X = vtx(:,1); Y = vtx(:,2); Z = vtx(:,3);
U = real(jmsh{1}); V = real(jmsh{2});  W = real(jmsh{3}); 
quiver3(X,Y,Z,U,V,W,'r')
hold on
plot(mesh,jnorm)
axis equal
title('|j|')
xlabel('X');   ylabel('Y');   zlabel('Z');
colorbar

% RESULTING TOTAL ELECTRIC FIELD
% Use Stratton-Chu representation

radiat = mshSquare(5e3,[4 4]);

% Esc = mathcal{I}j - matchal{K}m
% where j=-nxHtot is the solution above

I1 = 1i * k * 1/(4*pi)*integral(radiat.vtx, Gamma, Gxy,rest(v,1))... 
        + 1i/k * 1/(4*pi) * integral(radiat.vtx, Gamma,gradxGxy{1}, div(v));

I2 = 1i * k * 1/(4*pi)*integral(radiat.vtx, Gamma, Gxy,rest(v,2))... 
        + 1i/k * 1/(4*pi) * integral(radiat.vtx, Gamma, gradxGxy{2}, div(v));

I3 = 1i * k * 1/(4*pi)*integral(radiat.vtx, Gamma, Gxy,rest(v,3))... 
        + 1i/k * 1/(4*pi) * integral(radiat.vtx, Gamma, gradxGxy{3}, div(v));
    
Esc1 = I1*j;
Esc2 = I2*j;
Esc3 = I3*j;

Einc1 = PWE{1}(radiat.vtx);
Einc2 = PWE{2}(radiat.vtx);
Einc3 = PWE{3}(radiat.vtx);

Etot1 = Einc1 + Esc1;
Etot2 = Einc2 + Esc2;
Etot3 = Einc3 + Esc3;

EtotNorm = sqrt(real(Etot1).^2 + real(Etot2).^2 + real(Etot3).^2);


% Graphical representation
figure

vtx = radiat.vtx;
X = vtx(:,1); Y = vtx(:,2); Z = vtx(:,3);
U = real(Etot1); V = real(Etot2);  W = real(Etot3); 
quiver3(X,Y,Z,U,V,W,'r')
hold on
plot(radiat,EtotNorm)
axis equal
title('|j|')
xlabel('X');   ylabel('Y');   zlabel('Z');
colorbar
caxis([0,2])
%% Checking Stratton Chu for the incident field (Must give 0 outside)

% Kb = curl int_Gamma G(x,y)b(y)dy =  int_Gamma (Nabla_x G(x,y)) x b(y)dy
% x notin Gamma

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

j = nxHinc;
m = nxEinc;

E1test = I1*j - K1*m;
E2test = I2*j - K2*m;
E3test = I3*j - K3*m;

EtestNorm = sqrt(real(E1test).^2 + real(E2test).^2 + real(E3test).^2);
figure

vtx = radiat.vtx;
X = vtx(:,1); Y = vtx(:,2); Z = vtx(:,3);
U = real(E1test); V = real(E2test);  W = real(E3test); 
quiver3(X,Y,Z,U,V,W,'r')
hold on
radiat.vtx(:,3) = EtestNorm;
plot(radiat,EtestNorm)
axis equal
title('|j|')
xlabel('X');   ylabel('Y');   zlabel('Z');
colorbar

