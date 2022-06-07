%+========================================================================+
% Maxwell electric wave equation direct first-kind BIE
% (Extension by zero approach in the interior)
%+========================================================================+

% Cleaning
clear all
close all
clc

% Gypsilab path
addpath(genpath("../../"));

% Parameters
N   = 1e3; 
tol = 1e-3; 
gss = 3; 

% Spherical mesh
mesh = mshSphere(N,1);
Gamma  = dom(mesh,3);
figure
plot(mesh)
hold on
plotNrm(mesh)

% Finite elements
v = fem(mesh,'RWG');

% Setting a frequency suitable for testing
stp = mesh.stp;
k   = 1/stp(2);
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

disp('Assembling M')
tic
M = - 1/(4*pi) * integral(Gamma, Gamma,nx(v), gradxGxy, v, tol); 
toc

disp('Regularizing M')
tic
% + because grady
M = M + 1/(4*pi) * regularize(Gamma, Gamma,nx(v),'grady[1/r]', v);
toc

Id = integral(Gamma,v,v);

% Incident currents 
nxEinc = - integral(Gamma,v,v)\(integral(Gamma,nx(v),PWE)); 
nxHin = - integral(Gamma,v,v)\(integral(Gamma,nx(v),PWH)); 


disp('Solving integral equation')
j = -1i*k*(T \ (-1/2*Id + M)*nxEinc);


%%% SURFACE CURRENT
disp('~~~~~~~~~~~~~ Surface current ~~~~~~~~~~~~~')

% Mesh Interpolation
jmsh = feval(v, j,mesh);
jnorm    = sqrt(sum(abs(cell2mat(jmsh)).^2,2));

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

radiat = mshSquare(1e3,[4 4]);

% Esc = mathcal{I}j - matchal{K}m
% where j=-nxHsc is the solution above and m = nxEinc

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
    
Esc1 = I1*j - K1*nxEinc;
Esc2 = I2*j - K2*nxEinc;
Esc3 = I3*j - K3*nxEinc;

% Graphical representation
figure

EscNorm = sqrt(abs(Esc1).^2 + abs(Esc2).^2 + abs(Esc3).^2);

vtx = radiat.vtx;
X = vtx(:,1); Y = vtx(:,2); Z = vtx(:,3);
U = real(Esc1); V = real(Esc2);  W = real(Esc3); 
quiver3(X,Y,Z,U,V,W,'r')
hold on
plot(radiat,EscNorm)
axis equal
title('|j|')
xlabel('X');   ylabel('Y');   zlabel('Z');
colorbar



