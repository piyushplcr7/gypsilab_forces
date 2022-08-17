addpath('../openMsh');
addpath('../openDom');
addpath('../openFem');
addpath('../openHmx');
addpath('../openBmm');
addpath('../openFfm');
addpath('../miscellaneous');  
addpath('../Source_SS_quad');

% clc
disp("----------------------------------------- ")
disp("----------------------------------------- ")


tic

mesh = mshCube(600, [2 2 2]);
kappa = 2;

tol = 1e-6;


Omega = dom(mesh, 4);

m_Gamma = bnd(mesh)
Gamma = dom(m_Gamma, 3);

Vh = fem(m_Gamma, 'P0');


elt2dof = Gamma.msh.elt;

Ndof = size(Gamma.msh.vtx, 1);
Nelt = size(elt2dof, 1);

Adj = sparse((1:Nelt)', elt2dof(:, 1), 1, Nelt, Ndof) + ...
    sparse((1:Nelt)', elt2dof(:, 2), 1, Nelt, Ndof) + ...
    sparse((1:Nelt)', elt2dof(:, 3), 1, Nelt, Ndof);


[I, J, Case] = find(Adj * Adj');

corr = ~~(Adj * Adj');

corr = ~~kron(corr, ones(Gamma.gss));

toc

Gxy = @(X, Y) femGreenKernel(X, Y, '[exp(ikr)/r]', kappa);

disp(' ')
disp('Galerkin Matrix - Far field interactions - H-Matrix')
A0 = 1 / (4*pi) * integral(Gamma, Gamma, Vh, Gxy, Vh, tol, corr);

toc


disp(' ')
disp('Galerkin Matrix - Far field interactions - Full-Matrix')
A0full = 1 / (4*pi) * integral(Gamma, Gamma, Vh, Gxy, Vh, corr);

toc


disp(' ')
disp('Galerkin Matrix - Regularized')
A1 = 1 / (4*pi) * integral(Gamma, Gamma, Vh, Gxy, Vh) + ...
     1 / (4*pi) * regularize(Gamma, Gamma, Vh, '[1/r]', Vh);

toc


disp(' ')
disp('Galerkin Matrix - Sauter-Schwab Near Field Matrix')
kernel = @(x,y,z) exp(1i*kappa*sqrt(sum(z.^2 ,2))) .* sqrt(1./ sum(z.^2 ,2) ) /4./pi;
A3 = panel_assembly(m_Gamma, kernel, Vh, Vh, I, J);

toc

% disp(' ')
% disp('Galerkin Matrix - Sauter-Schwab Full Matrix')
kernel = @(x,y,z) exp(1i*kappa*sqrt(sum(z.^2 ,2))) .* sqrt(1./ sum(z.^2 ,2) ) /4./pi;
A2 = panel_oriented_assembly(m_Gamma, kernel, Vh, Vh);
toc




A = A0 + A3;
Afull = A0full + A3;

norm(A1 - full(A0))

norm(A1 - A0full)

norm(A1 - A3 - full(A0))

norm(A1 - full(A))

norm(A1 - Afull)

norm(Afull - A2)
