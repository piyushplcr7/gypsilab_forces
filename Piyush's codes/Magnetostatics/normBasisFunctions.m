% Magnetostatics transmission problem with a source

addpath(genpath("../../"));
clear; clc; close all;

ivals = 5:12;
normRWG0 = 0*ivals;
normRWG=normRWG0;
normTEST = normRWG;
normCURL = normTEST;
hvals = normRWG0;


Nvals = size(ivals,2);

for i = 1:Nvals
N = 2^ivals(i);
%% SOLUTION DOMAIN
% Cube size and position
L = 2*[1 1 1];
T = [5 5 3];

% Cube domain
%bndmesh = bndmeshCubeTranslated(N,L,T);

% Spherical domain
%bndmesh = mshSphere(N,1);
%bndmesh = bndmesh.translate(T);

mesh = mshCube(N,L);
mesh = mesh.translate(T);
%mesh = mesh.sub(1);
bndmesh = mesh.bnd;

% Mesh size
hvals(i) = sqrt(mean(bndmesh.ndv,1));

%% BEM
Gamma = dom(bndmesh,3);

% BEM spaces

P1 = fem(bndmesh,'P1');
% Div conforming with div0 constraint -> Neumann trace
RWG0 = nxgrad(P1); 
% Div conforming space 
RWG = fem(bndmesh,'RWG');

%% Matrices containing norm^2

SL_RWG0_RWG0 = single_layer(Gamma,RWG0,RWG0); % H-1/2(divT,T) norm for DIV0 functions

% % H-1/2(divT,T) norm for RWG functions
SL_RWG_RWG = single_layer(Gamma,RWG,RWG); 
% For -1/2 norm of surface divergence
Norm2div = single_layer(Gamma,RWG.div,RWG.div); 

% Matrix for curl norm
SL_NED_NED = single_layer(Gamma,RWG.nx,RWG.nx);

normRWG0(i) = sqrt(max(diag(SL_RWG0_RWG0)))
normRWG(i) = sqrt(max(diag(SL_RWG_RWG + Norm2div)))
normTEST(i) = sqrt(max(diag(SL_RWG_RWG)))
normCURL(i) = sqrt(max(diag(SL_NED_NED)))

end

loglog(hvals,normRWG0,'-^');
mdl1 = fitlm(log(hvals),log(normRWG0));
m1 = mdl1.Coefficients.Estimate(2);
disp(m1);
c1 = mdl1.Coefficients.Estimate(1);
estimate1 = exp(m1 * log(hvals) + c1);
hold on;

loglog(hvals,estimate1,'-.')

cutoff = 1;
loglog(hvals,normRWG,'-*');
mdl2 = fitlm(log(hvals(cutoff:end)),log(normRWG(cutoff:end)));
m2 = mdl2.Coefficients.Estimate(2);
c2 = mdl2.Coefficients.Estimate(1);
disp(m2);
estimate2 = exp(m2 * log(hvals(cutoff:end)) + c2);

loglog(hvals(cutoff:end),estimate2,'--')

loglog(hvals,normTEST,'-+');
mdl3 = fitlm(log(hvals),log(normTEST));
m3 = mdl3.Coefficients.Estimate(2);
disp(m3);
c3 = mdl3.Coefficients.Estimate(1);
estimate3 = exp(m3 * log(hvals) + c3);
hold on;

loglog(hvals,estimate3,'-.')

loglog(hvals,normCURL,'-+');
mdl4 = fitlm(log(hvals),log(normCURL));
m4 = mdl4.Coefficients.Estimate(2);
disp(m4);
c4 = mdl4.Coefficients.Estimate(1);
estimate4 = exp(m4 * log(hvals) + c4);
hold on;

loglog(hvals,estimate4,'-.')