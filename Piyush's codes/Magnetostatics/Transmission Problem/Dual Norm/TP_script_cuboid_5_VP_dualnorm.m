% New transmission problem script
delete(gcp('nocreate'));
addpath(genpath("../../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 4;
mu0 = 2;
vals = 5:12;
Nvals = size(vals,2);
forces_mst = zeros(Nvals,3);
% indices for velocity fields go from 0 to kappa-1
kappa = 3;
shape_derivatives_mst = zeros(Nvals,3 * kappa^3);
shape_derivatives_bem = shape_derivatives_mst;

% Map of velocity field to index
abc_alpha = zeros(3 * kappa^3,4);

for a = 0:kappa-1
    for b = 0:kappa-1
        for c = 0:kappa-1
            for alpha = 0:2
                idx = a + kappa * b + kappa^2 * c + kappa^3 * alpha + 1;
                abc_alpha(idx,:) = [a b c alpha];
            end
        end
    end
end

testidx = 1 + kappa * 1 + kappa^2 * 1 + kappa^3 * 0 + 1;

% abc_alpha = [1 1 1 0];
% abc_alpha = repelem(abc_alpha,4,1);

% Number of fields
Nfields = size(abc_alpha,1);

forces_bem = forces_mst;
torques_mst = forces_mst;
torques_bem = forces_mst;
hvals = 0*vals;

for i = 1:Nvals
    N = 2^vals(i);
    disp(N);
    %% SOLUTION DOMAIN
    % Cube size and position
    L = [3 1 1];
    T = [2 1 3];
    
    % Cube domain
    %bndmesh = bndmeshCubeTranslated(N,L,T);
    
    % Spherical domain
%     bndmesh = mshSphere(N,1);
%     bndmesh = bndmesh.translate(T);
    
    mesh = mshCube(N,L);
    mesh = mesh.translate(T);
%     %mesh = mesh.sub(1);
    bndmesh = mesh.bnd;
    
    % Mesh size
    hvals(i) = sqrt(mean(bndmesh.ndv,1));
    
    Gamma = dom(bndmesh,3);
    normals = Gamma.qudNrm;
    
    %% Source
    N_src = N;
    R0 = 2;
    r0 = .5;
    [J,mesh_src] = get_torus_source(N_src,R0,r0);
    omega_src = dom(mesh_src,3);
    
    %% Solving the transmission problem
    % These are traces from the exterior
    % Psi lies in the space nxgradP1 and g lies in the space NED
    [Psi,g] = solveTransmissionProblem(bndmesh,J,omega_src,mu,mu0);
    % Interior traces
    Psi_in = mu/mu0 * Psi;
    g_in = g;
    
    %% Computing the MST based force and torque
    
    % Force computation
    NED = fem(bndmesh,'NED'); 
    P1 = fem(bndmesh,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0 = nxgrad(P1); 
    RWG = fem(bndmesh,'RWG');
    
    % Bn = curlA.n = curlTg
    Bn = reconstruct(g,Gamma,NED.curl);
    % Ht = nx(Hxn) = mu_e^-1 nxPsi
    Psivals = reconstruct(Psi,Gamma,DIV0);
    Ht = mu0^(-1) * cross(normals,Psivals,2);
    Ht = vecnorm(Ht,2,2);

    for fieldID = 1:Nfields
        a = abc_alpha(fieldID,1);
        b = abc_alpha(fieldID,2);
        c = abc_alpha(fieldID,3);
        alpha = abc_alpha(fieldID,4);
        [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);
        shape_derivatives_mst(i,fieldID) = ShapeDervTpVol(Gamma,Bn,Ht,mu0,mu,Vel);
    end
    

    %% Computing forces and torques using BEM shape derivative

    % Projecting Psi to RWG
    Psi_RWG = proj(Psivals,Gamma,RWG);

    [Vel1,DVel1] = getTransVelDVel([1 0 0]);
    [Vel2,DVel2] = getTransVelDVel([0 1 0]);
    [Vel3,DVel3] = getTransVelDVel([0 0 1]);

    nf1 = ForceSdBemTP(bndmesh,Psi,g,J,omega_src,Vel1);
    nf2 = ForceSdBemTP(bndmesh,Psi,g,J,omega_src,Vel2);
    nf3 = ForceSdBemTP(bndmesh,Psi,g,J,omega_src,Vel3);

    forces_bem(i,:) = [nf1 nf2 nf3]
    
    [Vel,DVel] = getCosVelDVel(1,1,1,1);

    shape_derivatives_bem(i,:) = SdBemTPVP_dualnorm(bndmesh,Psi_RWG,g,J,omega_src,mu0,mu,abc_alpha);

end
