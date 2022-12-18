% Solving curl curl problem

%\int_{\Omega} nu curlA.curlA' dx + \int_{\Omega} grad psi. A' dx +
%\int_{\Omega} grad phi. A dx = \int_{\Omega} J . A' dx

addpath(genpath("../../"));

clear; clc;

N = 50;
L = 10*[1 1 1];

mesh = mshCube(N,L);
bnd_mesh = mesh.bnd;

Vcurl = fem(mesh,'NED');
V = fem(mesh,'P1');

% Dirichlet boundary conditions
Vcurl0 = dirichlet(Vcurl,mesh.bnd);
V0 = dirichlet(V,mesh.bnd);

qud_order = 4;
omega = dom(mesh,qud_order);

% Reluctivity
nu = 1;

curlcurlmat = nu * integral(omega,Vcurl0.curl,Vcurl0.curl);
mixmat = integral(omega,Vcurl0,V0.grad);
N_Vcurl0 = size(curlcurlmat,1);
N_V0 = size(mixmat,2);
zeromat = zeros(N_V0,N_V0);

sysmat = [curlcurlmat mixmat; mixmat' zeromat];

% Test source current- solenoid around origin
% 
R0 = 2;
r0 = .5;

interior_torus = @(x,y,z) sqrt( (x-x*R0./(sqrt(x.^2+y.^2))).^2 + (y-y*R0./(sqrt(x.^2+y.^2))).^2 + z.^2 )<r0;
torus_tangent = @(x,y,z) interior_torus(x,y,z).*[-y x z*0]./(sqrt(x.^2+y.^2));

J0 = 1;

J = @(x,y,z) J0*torus_tangent(x,y,z); 

JJ = @(X) J(X(:,1),X(:,2),X(:,3));
%J = [1;0;0];

% Want to write an integral like integral(omega,Vcurl0,J)
%pre_rhs_top = integral(omega,JJ,Vcurl0);
% Performing manual integration for pre_rhs_top
[X,W] = omega.qud;
uqmat = Vcurl0.uqm(omega);

JX = J(X(:,1),X(:,2),X(:,3));

pre_rhs_top = sum(W.*(JX(:,1).*uqmat{1} + JX(:,2).*uqmat{2} + JX(:,3).*uqmat{3}),1)';


% Explicit computation possible only for specific field
N_temp = size(pre_rhs_top,1);
rhs_top = eye(N_temp,N_temp)*pre_rhs_top;
rhs_bot = zeros(N_V0,1);
% Not computing explicitly because of the simple field
rhs_vec = [rhs_top; rhs_bot];
%rhs_vec = 0*rhs_vec;
%rhs_vec(1) = 1;

sol = sysmat\rhs_vec;

sol_vecpot = sol(1:N_Vcurl0);
sol_psi = sol(N_Vcurl0+1:end);

% Projection to full fem spaces without Dirichlet BC
P_Curl_Curl0 = elimination(Vcurl,mesh.bnd);
sol_vecpot_full = P_Curl_Curl0 * sol_vecpot;

P_V_V0 = elimination(V,mesh.bnd);
sol_psi_full = P_V_V0 * sol_psi;

% Getting the quadrature nodes and weights
[X,W,elt2qud] = omega.qud;

% uqmats contain the value of the basis function at the quadrature points
% (Rows denote quadrature points and columns denote the basis functions)
Vcurl_uqmat = Vcurl.uqm(omega);
V_uqmat = V.uqm(omega);

N_Vcurl = Vcurl.ndof;
N_V = V.ndof;

% Extracting the 3 components of the NED basis functions at the quadrature
% points
Id = eye(N_Vcurl,N_Vcurl);
Vcurl_uq1 = Vcurl_uqmat{1} * Id;
Vcurl_uq2 = Vcurl_uqmat{2} * Id;
Vcurl_uq3 = Vcurl_uqmat{3} * Id;

% We need to multiply by the coefficients to get the solution values at
% quadrature points
sol_vecpot_uqmat1 = bsxfun(@times,sol_vecpot_full',Vcurl_uq1);
sol_vecpot_uqmat2 = bsxfun(@times,sol_vecpot_full',Vcurl_uq2);
sol_vecpot_uqmat3 = bsxfun(@times,sol_vecpot_full',Vcurl_uq3);

sol_vecpot_vec1 = sum(sol_vecpot_uqmat1,2);
sol_vecpot_vec2 = sum(sol_vecpot_uqmat2,2);
sol_vecpot_vec3 = sum(sol_vecpot_uqmat3,2);

% Making a vector for visualization. 
%Nqud = size(sol_vecpot_uqmat3,1);
%Nbas = size(sol_vecpot_uqmat3,2);
%newsize = [Nqud*Nbas,1];
%sol_vecpot_vec1 = reshape(sol_vecpot_uqmat1,newsize);
%sol_vecpot_vec2 = reshape(sol_vecpot_uqmat2,newsize);
%sol_vecpot_vec3 = reshape(sol_vecpot_uqmat3,newsize);

%sol_vecpot_vec = [sol_vecpot_vec1 sol_vecpot_vec2 sol_vecpot_vec3];
%sol_vecpot_vec = vecnorm(sol_vecpot_vec,2,2);

% Need a vector of quadrature points, repeated Nbas no. of times
%Xrep = repmat(X,Nbas,1);

% Visualizing in a quiver3 plot
figure;
m = mesh.bnd;
H = patch('Faces',m.elt,'Vertices',m.vtx);
set(H,'FaceVertexCData',0,'FaceColor','flat');
set(H,'FaceAlpha', 0.2);
hold on;
quiver3(X(:,1),X(:,2),X(:,3),sol_vecpot_vec1,sol_vecpot_vec2,sol_vecpot_vec3);

% Visualizing the curl of the vector potential or the magnetic induction
curlVcurl = Vcurl.curl;

curlVcurl_uqmat = curlVcurl.uqm(omega);

N_curlVcurl = curlVcurl.ndof;

% Extracting the 3 components of the NED basis functions at the quadrature
% points
Id = eye(N_curlVcurl,N_curlVcurl);
curlVcurl_uq1 = curlVcurl_uqmat{1} * Id;
curlVcurl_uq2 = curlVcurl_uqmat{2} * Id;
curlVcurl_uq3 = curlVcurl_uqmat{3} * Id;

% We need to multiply by the coefficients to get the solution values at
% quadrature points
sol_magind_uqmat1 = bsxfun(@times,sol_vecpot_full',curlVcurl_uq1);
sol_magind_uqmat2 = bsxfun(@times,sol_vecpot_full',curlVcurl_uq2);
sol_magind_uqmat3 = bsxfun(@times,sol_vecpot_full',curlVcurl_uq3);

sol_magind_vec1 = sum(sol_magind_uqmat1,2);
sol_magind_vec2 = sum(sol_magind_uqmat2,2);
sol_magind_vec3 = sum(sol_magind_uqmat3,2);

figure;
m = mesh.bnd;
H = patch('Faces',m.elt,'Vertices',m.vtx);
set(H,'FaceVertexCData',0,'FaceColor','flat');
set(H,'FaceAlpha', 0.2);
hold on;
quiver3(X(:,1),X(:,2),X(:,3),sol_magind_vec1,sol_magind_vec2,sol_magind_vec3);

%quiver3(Xrep(:,1),Xrep(:,2),Xrep(:,3),sol_vecpot_vec(:,1),sol_vecpot_vec(:,2),sol_vecpot_vec(:,3));

% Visualizing the vector potential
%plot(mesh,sol_vecpot_full)

% Visualizing the magnetic field as the curl of the vector potential
% Need curl of the basis functions

