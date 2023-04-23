% Testing Calderon Identity for Laplace

addpath(genpath("../../"));
clear; clc; close all;

ivals = 5:12;
err1 = 0*ivals;
err2=err1;
hvals = err1;

Nvals = size(ivals,2);

for i = 1:Nvals

N = 2^ivals(i);

L = 2*[1 1 1];
T = [5 5 3];

% Cube domain
bndmesh = bndmeshCubeTranslated(N,L,T);

%bndmesh = mshSphere(N,1);
%bndmesh = bndmesh.translate([2 2 2]);
hvals(i) = sqrt(mean(bndmesh.ndv,1));

Gamma = dom(bndmesh,3);
P1 = fem(bndmesh,'P1');
P0 = fem(bndmesh,'P0');

% Defining the kernel for SL
Gxy = @(X,Y)femGreenKernel(X,Y,'[log(r)]',0); % 0 wave number
% Defining the kernel for DL
GradGn = cell(3,1);
GradGn{1} = @(X,Y)femGreenKernel(X,Y,'grady[log(r)]1',0);
GradGn{2} = @(X,Y)femGreenKernel(X,Y,'grady[log(r)]2',0);
GradGn{3} = @(X,Y)femGreenKernel(X,Y,'grady[log(r)]3',0);

% Evaluating the bilinear form for the Single Layer BIO
V = single_layer(Gamma,P0,P0);

% Evaluating the bilinear form for the Double Layer BIO
K = double_layer_laplace(Gamma,P0,P1);

% Hypersingular
W = single_layer(Gamma,nxgrad(P1),nxgrad(P1));

[X,Wt] = Gamma.qud;


%%
% Dirichlet trace
TD = 1./(4*pi * vecnorm(X,2,2));
%TD = sum(X,2);
%TD = X(:,1).^2 + X(:,2).^2 - 2*X(:,3).^2;

% Neumann trace
normals = Gamma.qudNrm;
normX3 = vecnorm(X,2,2).^3;
xdotnrm = dot(X,normals,2);
TN = -1/(4*pi) * xdotnrm./normX3;
%TN = sum(normals,2);
%TN = dot(2*[X(:,1) X(:,2) -2*X(:,3)],normals,2);

TD_P1_coeffs = proj(TD,Gamma,P1);
TN_P0_coeffs = proj(TN,Gamma,P0);

%%

M = mass_matrix(Gamma,P0,P1);

traces = [TN_P0_coeffs; TD_P1_coeffs];

vec1 = [V -0.5*M-K] * traces;

vec2 = [0.5*M'-K' -W] * traces;

disp(N)
err1(i) = max(abs(vec1))
err2(i) = max(abs(vec2))

end

loglog(hvals,err1,'-^');
mdl1 = fitlm(log(hvals),log(err1));
m1 = mdl1.Coefficients.Estimate(2);
disp(m1);
c1 = mdl1.Coefficients.Estimate(1);
estimate1 = exp(m1 * log(hvals) + c1);
hold on;

loglog(hvals,estimate1,'-.')

cutoff = 3;
loglog(hvals,err2,'-*');
mdl2 = fitlm(log(hvals(cutoff:end)),log(err2(cutoff:end)));
m2 = mdl2.Coefficients.Estimate(2);
c2 = mdl2.Coefficients.Estimate(1);
disp(m2);
estimate2 = exp(m2 * log(hvals(cutoff:end)) + c2);

loglog(hvals(cutoff:end),estimate2,'--')



