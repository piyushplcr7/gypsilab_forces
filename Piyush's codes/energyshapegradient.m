% inputs: mesh needs to be a boundary mesh as we need to compute the full
% gradients at the boundary

function sg = energyshapegradient(mesh,g,nu)

bnd_mesh = mesh.bnd;

% getting the state solution
psi = Dirichlet_DFK(bnd_mesh,g);
% getting the adjoint solution
p = bem_sg_adjoint(bnd_mesh,g);

% Creating the necessary BEM spaces
S1_Gamma = fem(bnd_mesh,'P1');
S0_Gamma = fem(bnd_mesh,'P0');
S1_Omega = fem(mesh,'P1');
S0_Omega = fem(mesh,'P0');

% Creating the integration domains
Gamma = dom(mesh,3);
Omega = dom(mesh,3);

% evaluating g at volume dofs
g_N = g(S1_Omega.dof);

% Evaluating the coefficients for gradient of g, which lies in the space P0
gradg = cell(3,1);
M00 = integral(Omega,S0_Omega,S0_Omega);
gradg{1} = M00\(integral(Omega,S0_Omega,grad(S1_Omega,1)) * g_N);
gradg{2} = M00\(integral(Omega,S0_Omega,grad(S1_Omega,2)) * g_N);
gradg{3} = M00\(integral(Omega,S0_Omega,grad(S1_Omega,3)) * g_N);

% Getting the P0_Omega to P0_Gamma trace operator
tr_opr = P0_vol_to_bnd(mesh);

% Getting gradg at the boundary dofs
gradgb = cell(3,1);
gradgb{1} = tr_opr * gradg{1};
gradgb{2} = tr_opr * gradg{2};
gradgb{3} = tr_opr * gradg{3};

% Getting coefficients for gradg.nu \in S0_Gamma
gradgvec = [gradgb{1} gradgb{2} gradgb{3}];
nuvec = nu(S0_Gamma.dof);
gradg_nu = sum(gradgvec.*nuvec,2);

% Getting coefficients for div nu \in S0_Gamma
nuvec = nu(S1_Omega.dof);
nu_1_vec = nuvec(:,1);
nu_2_vec = nuvec(:,2);
nu_3_vec = nuvec(:,3);
dnu1_x = M00\(integral(Omega,S0_Omega,grad(S1_Omega,1)) * nu_1_vec);
dnu2_y = M00\(integral(Omega,S0_Omega,grad(S1_Omega,2)) * nu_2_vec);
dnu3_z = M00\(integral(Omega,S0_Omega,grad(S1_Omega,3)) * nu_3_vec);
divnu = dnu1_x + dnu2_y + dnu3_z;
% Getting coefficients for div(g nu) = gradg.nu + g div nu
div_g_nu = gradg_nu + g(S0_Omega.dof).*divnu;
divgnu_bnd = tr_opr * div_g_nu;

% Calculating the individual terms for the shape gradient

% T1 
M00_bnd = integral(Gamma,S0_Gamma,S0_Gamma);
T1 = 0.5 * psi' * M00_bnd * gradg_nu;

% T6
T6 = -0.5 * p' * M00_bnd * gradg_nu;

% T5
% Defining the kernel for DL
GradGn = cell(3,1);
GradGn{1} = @(X,Y)femGreenKernel(X,Y,'grady[log(r)]1',0);
GradGn{2} = @(X,Y)femGreenKernel(X,Y,'grady[log(r)]2',0);
GradGn{3} = @(X,Y)femGreenKernel(X,Y,'grady[log(r)]3',0);
% Evaluating the bilinear form for the Double Layer BIO
K = -1/(2*pi)*integral(Gamma,Gamma,S0_Gamma,GradGn,ntimes(S1_Gamma));
% Regularization
K = K -1/(2*pi)*regularize(Gamma,Gamma,S0_Gamma,'grady[log(r)]',ntimes(S1_Gamma));
T5 = -p' * K * divgnu_bnd;

% T2
t1 = integral(Gamma,Gamma,S1_Gamma,GradGn{1},S1_Gamma);
t2 = integral(Gamma,Gamma,S1_Gamma,GradGn{2},S1_Gamma);
t3 = integral(Gamma,Gamma,S1_Gamma,GradGn{3},S1_Gamma);
vels = nu(S1_Gamma.dof);
NU1 = diag(vels(:,1));
NU2 = diag(vels(:,2));
NU3 = diag(vels(:,3));

t1 = t1 - regularize(Gamma,Gamma,S1_Gamma,'grady[1/r]',S1_Gamma);
t2 = t2 - regularize(Gamma,Gamma,S1_Gamma,'grady[1/r]',S1_Gamma);
t3 = t3 - regularize(Gamma,Gamma,S1_Gamma,'grady[1/r]',S1_Gamma);

t1 * NU1 + t2 * NU2 + t3 * NU3;


% Summing up terms to get the shape gradient
sg = T1 + T6 + T5;
end

function y = nu(X)
    y = X;
end
