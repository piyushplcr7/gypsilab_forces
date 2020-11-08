%% Liste des choses à faire : 
% - Maillage d'un tore -> done
% - Comprendre les éléments H(curl) conformes -> Nédélec
% - Regularization -> done


%% Hodge-Laplacian, T-formulation. 
% Validation with traces of harmonic functions. 
clear all
close all
clc

N = 500;
mesh = mshSphere(N,1);

Gamma = dom(mesh,3);

Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0);
dG{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
dG{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
dG{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);

Vh = fem(mesh,'RWG'); % Vectors Rao-Wilton-Galtson <=> RT0
Wh = fem(mesh,'P0'); 

IVh = integral(Gamma,Vh,Vh); % L^2 Mass RT0
IWh = integral(Gamma,Wh,Wh); % L^2 Mass P0

% Hodge Laplacian : block matrix of bilinear form
A = 1/(4*pi)*(integral(Gamma,Gamma,Vh,Gxy,Vh) + ...
    regularize(Gamma,Gamma,Vh,'[1/r]',Vh));
B = 1/(4*pi)*(integral(Gamma,Gamma,div(Vh),Gxy,Wh)+ ...
    regularize(Gamma,Gamma,div(Vh),'[1/r]',Wh));
BT = transpose(B);
Z = sparse([],[],[],size(B,2),size(B,2));
system = [[(A+A')/2, B]; [BT, Z]];



% Some harmonic vector field:
[X,Y,Z] = FunR3.XYZ;
u1 = X^2 - Y^2; 
u2 = Y^2 - Z^2;
u3 = Z^2 - X^2;
u = cell(1,3);
u{1} = u1; u{2} = u2; u{3} = u3;

divu = 2*X + 2*Y + 2*Z;
nxcurlu{1} = 2*Z;
nxcurlu{2} = 2*X;
nxcurlu{3} = 2*Y;

% Let (q,v) = T_D(u), i.e.
% q = Y_t(u) = nx(uxn)
% v = Y_D(u) = div(u)|Gamma
% We replace q and v by approximations of them in Vh and Wh respectively

% We use the formula (nx(uxn),v) = (nx(vxn),u)
% and nx(vxn) = v - (n,v)n
% (nxv,nxu)

Ytu = integral(Gamma,nxxn(Vh),u); % = int_{Gamma} mu(x) nx(uxn)(x) dx
q = IVh\Ytu;
YDu = integral(Gamma,Wh,divu);
v = IWh\YDu;

qPlot = feval(Vh,q,mesh);
figure
subplot(1,3,1)
plotOn(mesh,qPlot{1})
title('First component')
axis equal
subplot(1,3,2)
plotOn(mesh,qPlot{2})
title('Second component')
axis equal
subplot(1,3,3)
plotOn(mesh,qPlot{3})
title('Third component')
axis equal

vPlot = feval(Wh,v,mesh);
figure
plotOn(mesh,vPlot)
axis equal
title('div(u)|\Gamma')
% L1 : Vh -> R defined by L1(mu) = L11(mu) + L12(mu) + L13(mu)
% L11(mu) = 1/2(q,mu)

L11 = 1/2*Ytu;

% L12(mu) = -({Yt} curl_x Psi (qxn) , mu)
%         = iint [grad_x (G(x-y)), mu(x), qxn(y)] dx dy where [u,v,w] =
%         (uxv).w (cf pdf note)
%        
%         = - iint [mu, grad_y(G(x-y)), nxq(y)] dx dy
% (- sign for the 3 sign changes). 
nxKxn = -1/(4*pi)*(integral(Gamma,Gamma,Vh,dG,nx(Vh))...
    + regularize(Gamma,Gamma,Vh,'grady[1/r]',nx(Vh)));

L12 = -nxKxn*q;


% L13(mu) = ({Yt} Gamma(v),mu)
%         =  iint (G(x-y) v(y) n(y)).mu(x) 
Vn = 1/(4*pi)*(integral(Gamma,Gamma,Vh,Gxy,ntimes(Wh)) +...
    regularize(Gamma,Gamma,Vh,'[1/r]',ntimes(Wh)));

L13 = Vn*v;
% L1 : 
L1 = L11 + L12 + L13;


% L2 : Vh -> R defied by L2(beta) = L21(beta) + L22(beta)
% L21(beta) = 1/2(v,beta)
L21 = 1/2*YDu;

% L22(beta) = -(div(Gamma(v)),beta)
%           = - iint grad_x(G(x-y)).n(y)v(y) beta(x)
%           =   iint grad_y(G(x-y)).n(y)v(y) beta(x)

K = 1/(4*pi)*(integral(Gamma,Gamma,Wh,dG,ntimes(Wh))...
    + regularize(Gamma,Gamma,Wh,'grady[1/r]',ntimes(Wh)));
% Kij   =      phi_i(x) d_x G(x-y) n(y) psi_j(y)
L22 = K*v; 
L2 = L21 + L22;

L = [L1;L2];
U = gmres(system,L);


f = U(1:size(Vh));
alpha = U((size(Vh)+1):end);

alphaPlot = feval(Wh,alpha,mesh);
figure;
plotOn(mesh,alphaPlot);

% Theoretical result :
% if alpha = u|Gamma . n

alphaTheo = IWh\integral(Gamma,ntimes(Wh),u);
alphaTheoPlot = feval(Wh,alphaTheo,mesh);
figure;
plotOn(mesh,alphaTheoPlot);

% f = nxcurl(u)

fTheo = IVh\integral(Gamma,Vh,nxcurlu);
plot(fTheo - f);


plotOn(mesh,alphaTheoPlot - alphaPlot);

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


