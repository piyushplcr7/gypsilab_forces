close all;
clear all;


%% Geometry
% Create a non radially symetric mesh mesh whose boundary has two connected components:
N = 30; rad1 = 1; rad2 = 2;
m_Omega = mshAssymetricRing(N,rad1,rad2);
m_Gamma = m_Omega.bnd; % boundary of the mesh
c = m_Gamma.ctr;
ind = c(:,1).^2 + c(:,2).^2 > (rad1 + rad2)^2/4;
m_Gamma1 = m_Gamma.sub(ind);
m_Gamma2 = setdiff(m_Gamma,m_Gamma1);

figure;
plot(m_Omega);
hold on;
plot(m_Gamma1,'g');
plot(m_Gamma2,'r');
axis off;
title('Mesh and components of the boundary');

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
figure
graph(S1_Omega,Gh)
view(45,10)
axis equal
axis off
title('Function Gh')

%% Computation of the electrostatic potential

% Linear form
rhs = -integral(Omega,grad(S10_Omega),grad(S1_Omega))*Gh;
% Bilinear form:
a = integral(Omega,grad(S10_Omega),grad(S10_Omega));
% Solution:
wh = a\rhs;

P = elimination(S1_Omega,m_Gamma);
Uh = P*wh + Gh;

figure
graph(S1_Omega,Uh);
view(45,10)
axis equal
axis off
title('Electrostatic potential')
axis tight;

%% Computation of the Electric field


S0_Omega = fem(m_Omega,'P0');
M0 = integral(Omega,S0_Omega,S0_Omega);
Omega = dom(m_Omega,gss); % We need just one Gauss point since everything is constant
E0 = cell(3,1);
E0{1} = M0\(integral(Omega,S0_Omega,grad(S1_Omega,1))*Uh);
% int_{Omega} psi_h \partial_x(phi_h)
E0{2} = M0\(integral(Omega,S0_Omega,grad(S1_Omega,2))*Uh);
E0{3} = M0\(integral(Omega,S0_Omega,grad(S1_Omega,3))*Uh);

figure
title('Electric field')
plot(m_Omega,'w');
hold on
c = m_Omega.ctr;
quiver(c(:,1),c(:,2),E0{1},E0{2},'r')

%% Computation of the normal trace


S0_Gamma2 = fem(m_Gamma2,'P0');
gss = 1;


% Compute the gradient of uh at the center of each element, and restrict to
% Gamma2

[tr_Gamma_P0,S0_Gamma] = trace_P0P1(S0_Omega);
P = restriction(S0_Gamma,m_Gamma2);
tr_Gamma2_P0 = P*tr_Gamma_P0;


gradUh1 = tr_Gamma2_P0*(E0{1});
gradUh2 = tr_Gamma2_P0*(E0{2});
gradUh3 = tr_Gamma2_P0*(E0{3});

% Scalar product by normals:
Nrm = m_Gamma2.nrm;
dnU = Nrm(:,1).*gradUh1 + Nrm(:,2).*gradUh2 + Nrm(:,3).*gradUh3;

figure
plot(m_Omega,'w');
hold on
surf(S0_Gamma2,dnU)

%% Net force by engineering formula:

Gamma2 = dom(m_Gamma2,1);
I = integral(Gamma2,S0_Gamma2,ntimes(S0_Gamma2));
Fengineers = cell(3,1);
for i = 1:3
    Fengineers{i} = 1/2*dnU'*I{i}*dnU;
end
disp(Fengineers{1}); % By symmetry, F_y = 0. 

%% Net force by volume formula:

E2 = E0{1}.^2 + E0{2}.^2 + E0{3}.^2;
E2xE = cell(3,1);
Fvol = cell(3,1);
for i = 1:3
    E2xE{i} = E2.*E0{i};
    Fvol{i} = 1/2*sum(m_Omega.ndv.*E2xE{i});
end



