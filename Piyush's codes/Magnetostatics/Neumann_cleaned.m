% curl curl Neumann problem
addpath(genpath("../../"));
clear; clc; close all;

ivals  = 6:10;
Nvals = size(ivals,2);
Hmherrs = ivals*0;
hvals = Hmherrs;

for i = 1:Nvals 
N = 2^ivals(i);
%% SOLUTION DOMAIN

% Cube domain
%bndmesh = bndmeshCubeTranslated(N,2*[1 1 1],[5 5 3]);

% Spherical domain
bndmesh = mshSphere(N,1);
bndmesh = bndmesh.translate([5 5 3]);

% Mesh size
%hvals(i) = sqrt(mean(bndmesh.ndv,1));

% BEM spaces
NED = fem(bndmesh,'NED'); % curl conforming -> Dirichlet trace
P1 = fem(bndmesh,'P1');
RWG0 = nxgrad(P1); % div conforming, div0 -> Neumann trace
RWG = fem(bndmesh,'RWG');
curlNED = NED.curl;

order = 7;

Gamma = dom(bndmesh,order);
hvals(i) = sqrt(mean(bndmesh.ndv,1));

%% Generating traces 
[TDA,TNA,curlAdotn] = getTraces1(Gamma);
%[TDA,TNA,curlAdotn] = getTraces2(Gamma);
% [TDA,TNA,curlAdotn,mesh_src] = getTracesTorusSource(Gamma); 

%% Solving 
% Solving for Dirichlet Trace using the second equation
TDA_sol_coeffs = Neumann_DFK(TNA,NED,'div0',order);
% TDA_sol_coeffs = Neumann_DFK_discrete(TNA,NED,'div0',order);

% Solving for Neumann Trace using the first equation
TNA_sol_coeffs = Dirichlet_DFK(TDA,RWG0,order);

% Solving for the Neumann Trace using the second equation
TNA_sol_coeffs_indirect = Dirichlet_DSK(TDA,RWG,order);

%% Computing the error for Neumann Problem
err_coeff = proj(TDA,Gamma,NED)-TDA_sol_coeffs;
disp(N)
% MinusHalfErrCurlg = err_coeff' * single_layer(Gamma,RWG.div,RWG.div) * err_coeff
Hmherrs(i) = err_coeff' * single_layer(Gamma,RWG.div,RWG.div) * err_coeff

% 
l2err = err_coeff' * mass_matrix(Gamma,NED,NED) * err_coeff

curlTg = reconstruct(TDA_sol_coeffs,Gamma,curlNED);
%norm(curlAdotn - curlTg)/norm(curlAdotn)
%norm(curlAdotn + curlTg)/norm(curlAdotn)


%% Visualizing the solution
if true
%     cutoff = 50;
%     %plot((curlAdotn + curlTg)./curlAdotn)
%     plot(curlAdotn(1:cutoff));
%     hold on;
%     plot(-curlTg(1:cutoff),'red');
% 
    % Plot without cutoff
    plot(curlAdotn);
    hold on;
    plot(-curlTg,'red');


    [X,~] = Gamma.qud;
    figure;
    quiver3wrapper(X,TDA,'blue');
    hold on;
    quiver3wrapper(X,reconstruct(TDA_sol_coeffs,Gamma,NED),'red');

    % Incorrect theoretically
    %quiver3wrapper(X,reconstruct(-TDA_sol_coeffs,Gamma,RWG),'red');
    title('Dirichlet trace computed with DFK');
    
    figure;
    quiver3wrapper(X,TNA,'blue');
    hold on;
    quiver3wrapper(X,reconstruct(TNA_sol_coeffs,Gamma,RWG0),'red');
    title('Neumann trace computed with DFK');

    figure;
    quiver3wrapper(X,TNA,'blue');
    hold on;
    quiver3wrapper(X,reconstruct(TNA_sol_coeffs_indirect,Gamma,RWG),'red');
    title('Neumann trace computed with DSK');
end
%loglog(hvals,Hmherrs);
end