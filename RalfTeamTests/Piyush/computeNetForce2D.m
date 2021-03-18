function [FBEM,Fvol,Feng,NvtxVol,NvtxBound,h] = computeNetForce2D(N)


%% Geometry
% Create a non radially symetric mesh mesh whose boundary has two connected components:

rad1 = 1; rad2 = 2;
m_Omega = mshAssymetricRing(N,rad1,rad2);
m_Gamma = m_Omega.bnd; % boundary of the mesh
c = m_Gamma.ctr;
indGamma1 = c(:,1).^2 + c(:,2).^2 > (rad1 + rad2)^2/4;
m_Gamma1 = m_Gamma.sub(indGamma1);
m_Gamma2 = setdiff(m_Gamma,m_Gamma1);


% Domains of integration
gss = 3; % Number of Gauss points
Omega = dom(m_Omega,gss);

%% Variational setting

% Finite element spaces
S1_Omega = fem(m_Omega,'P1'); % Piecewise linear functions on Omega. 
S10_Omega= dirichlet(S1_Omega,m_Gamma); % Pw linear Dirichlet condition

% Construct a function G on Omega which is equal to 1 on Gamma1 and 0 on
% Gamma2. 

P= restriction(S1_Omega,m_Gamma2); % operator restriction on Gamma2
Q = P'; % Adjoint operator: extension by 0. 
Gh= Q*ones(m_Gamma2.nvtx,1);


%% Computation of the electrostatic potential

% Linear form
rhs = -integral(Omega,grad(S10_Omega),grad(S1_Omega))*Gh;
% Bilinear form:
a = integral(Omega,grad(S10_Omega),grad(S10_Omega));
% Solution:
wh = a\rhs;

P = elimination(S1_Omega,m_Gamma);
Uh = P*wh + Gh;



%% Computation of the Electric field


S0_Omega = fem(m_Omega,'P0');
M0 = integral(Omega,S0_Omega,S0_Omega);
Omega = dom(m_Omega,gss); % We need just one Gauss point since everything is constant
E0 = cell(3,1);
E0{1} = M0\(integral(Omega,S0_Omega,grad(S1_Omega,1))*Uh);
E0{2} = M0\(integral(Omega,S0_Omega,grad(S1_Omega,2))*Uh);
E0{3} = M0\(integral(Omega,S0_Omega,grad(S1_Omega,3))*Uh);

%% Computation of the normal trace


gss = 1;


S0_Omega = fem(m_Omega,'P0');
M0 = integral(Omega,S0_Omega,S0_Omega);
Omega = dom(m_Omega,gss); % We need just one Gauss point since everything is constant

S0_Gamma2 = fem(m_Gamma2,'P0');


% Compute the gradient of uh at the center of each element, and restrict to
% Gamma2

[tr_Gamma_P0,S0_Gamma] = trace_P0P1(S0_Omega);
P = restriction(S0_Gamma,m_Gamma2);
tr_Gamma2_P0 = P*tr_Gamma_P0;


gradUh1 = tr_Gamma2_P0*(M0\integral(Omega,S0_Omega,grad(S1_Omega,1))*Uh);
gradUh2 = tr_Gamma2_P0*(M0\integral(Omega,S0_Omega,grad(S1_Omega,2))*Uh);
gradUh3 = tr_Gamma2_P0*(M0\integral(Omega,S0_Omega,grad(S1_Omega,3))*Uh);

% Scalar product by normals:
Nrm = m_Gamma2.nrm;
dnU = Nrm(:,1).*gradUh1 + Nrm(:,2).*gradUh2 + Nrm(:,3).*gradUh3;


%% Net force

Gamma2 = dom(m_Gamma2,1);
I = integral(Gamma2,S0_Gamma2,ntimes(S0_Gamma2));
Feng = 1/2*dnU'*I{1}*dnU;


E2 = E0{1}.^2 + E0{2}.^2 + E0{3}.^2;
Fvol = 1/2*sum(m_Omega.ndv.*E2.*E0{1});


NvtxVol= m_Omega.nvtx;

%% BEM approach

%% Geometry
% Create a non radially symetric mesh mesh whose boundary has two connected components:

% Domains of integration
gss = 3; % Number of Gauss points
Gamma = dom(m_Gamma,gss);

Gxy = @(X,Y)femGreenKernel(X,Y,'[log(r)]',0); % 0 wave number
S1_Gamma = fem(m_Gamma,'P0');
S1_Gamma2 = fem(m_Gamma2,'P0');

% Single layer potential
V = -1/(2*pi)*integral(Gamma,Gamma,S1_Gamma,Gxy,S1_Gamma);
V = V + (-1/(2*pi))*regularize(Gamma,Gamma,S1_Gamma,'[log(r)]',S1_Gamma);

B = integral(Gamma,S1_Gamma);
sys = [V B;B' 0];

% right hand side:
P = restriction(S1_Gamma,m_Gamma2);
F = integral(Gamma,S1_Gamma,S1_Gamma)*(P'*ones(size(P,1),1));
rhs = [F;0];
sol = sys\rhs;
lambda = sol(1:end-1);
dnu = P*lambda;

I = integral(Gamma2,S1_Gamma2,ntimes(S1_Gamma2));
FBEM = 1/2*dnu'*I{1}*dnu;

NvtxBound = m_Gamma.nvtx;

h = m_Omega.stp;
h = h(3);

end

