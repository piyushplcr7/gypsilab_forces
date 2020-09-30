%% Liste des choses à faire : 
% - Maillage d'un tore.
% - Comprendre les éléments H(curl) conformes
% - Regularization



%% Hodge Laplacian : T formulation


% Système à discrétiser : 
% a(mu,eta) = \int_{\Gamma x \Gamma} mu(x) G(x,y) eta(y)      mu, eta H(div)
% b(alpha,eta) = \int_{\Gamma x \Gamma} alpha G(x,y) div(eta)  alpha
% H^{-1/2}, eta H(div)

% H(div) RT_0
% H^{-1/2} P0 

% Stabilité inf sup. 
 

mesh = bnd(mshCube(100,[1 1 1]));
% plot(mesh);

Gamma = dom(mesh,3);

Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0);


Vh = fem(mesh,'RWG'); % Vectors Rao-Wilton-Galtson <=> RT0
Wh = fem(mesh,'P0');

A = 1/(4*pi)*integral(Gamma,Gamma,Vh,Gxy,Vh); % + regularization
B = 1/(4*pi)*integral(Gamma,Gamma,div(Vh),Gxy,Wh);
BT = transpose(B);

M11 = 1/(4*pi)*(integral(Gamma,Gamma,Vh,Gxy,Vh) + ...
    integral(Gamma,Gamma,div(Vh),Gxy,div(Vh))); 



M22 = 1/(4*pi)*integral(Gamma,Gamma,Wh,Gxy,Wh);
M = [M11 0*B; 0*BT M22];

Z = sparse([],[],[],size(B,2),size(B,2));


system = [[A, B]; [BT, Z]];

L = zeros(size(system,1),1);

U = system\L;

[P,D] = eig(full(M\system));
d = diag(D);

figure;
semilogy(sort(abs(d)),'*');

% We find the same result as in Claeys and Hiptmair



%% Hodge-Laplacian : N formulation.

clear all
close all
clc;

mesh = bnd(mshCube(100,[1 1 1]));
% plot(mesh);

Gamma = dom(mesh,3);
Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0);

Xh = fem(mesh,'NED');

A = integral(Gamma,Gamma,curl(Xh),Gxy,curl(Xh));

% curl_Gamma(p) = div(-n x p)




%% Dirac operator


clear all
close all
clc;

mesh = bnd(mshCube(100,[1 1 1]));
% plot(mesh);

Gamma = dom(mesh,3);
Vh = fem(mesh,'RWG');

