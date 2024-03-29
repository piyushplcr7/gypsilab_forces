% Script for constant Magnetization

addpath(genpath("../../../../"));
clear; clc; close all;
format long;

mu = 1;
vals = 5:13;
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
    % bndmesh = mshSphere(N,1);
    % bndmesh = bndmesh.translate(T);
    
    mesh = mshCube(N,L);
    mesh = mesh.translate(T);
%     %mesh = mesh.sub(1);
    bndmesh = mesh.bnd;
    
    % Mesh size
    hvals(i) = sqrt(mean(bndmesh.ndv,1));
    
    Gamma = dom(bndmesh,3);
    normals = Gamma.qudNrm;

    % Function to determine the interior of the magnet
    %interior = @(X) (vecnorm(X-T,2,2) < 1);
    
    %% Source
    N_src = N;
    R0 = 2;
    r0 = .5;
    [J_orig,mesh_src] = get_torus_source(N_src,R0,r0);
    omega_src = dom(mesh_src,3);
    
    %% Solving the Magnet problem
    % Constant Magnetization function
    M = @(X) ones(size(X,1),1) * [1 0 0];
    % Modifying J source
    J = @(x,y,z) 1 * J_orig(x,y,z);

    % These are traces from the exterior
    % Psi lies in the space nxgradP1 and g lies in the space NED
    [TnAJ,TdAJ,TnAM,TdAM] = solveMagnetProblemSimplified(Gamma,J,omega_src,mu,mu,M);

    % Projecting TnAM from nxgradP1 to RWG
    P1 = fem(bndmesh,'P1');
    RWG = fem(bndmesh,'RWG');
    PsivalsM = reconstruct(TnAM,Gamma,nxgrad(P1));
    TnAM_RWG = proj(PsivalsM,Gamma,RWG);

    PsivalsJ = reconstruct(TnAJ,Gamma,nxgrad(P1));
    TnAJ_RWG = proj(PsivalsJ,Gamma,RWG);

    %% Computing the B field for verification
%     plot_field_magnet(TdAM,TnAM_RWG,bndmesh,J,omega_src,mu0,interior);
%     figure;
%     plot_field_magnet(TdAJ,TnAJ_RWG,bndmesh,J,omega_src,mu0,interior);

    %% Computing forces and torques from Equivalent current model

    % Force computation
    NED = fem(bndmesh,'NED'); 
    P1 = fem(bndmesh,'P1');
    P0 = fem(bndmesh,'P0');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0 = nxgrad(P1); 
    RWG = fem(bndmesh,'RWG');
    
    % Testing the condition number
%     testmat = integral(Gamma,DIV0,DIV0);
%     rcond(full(testmat))

    [X,W] = Gamma.qud;

    Mvals = M(X);
    Mxn = cross(Mvals,normals,2);

    Mdotn = dot(Mvals,normals,2);
    Mdotn_coeffs = proj(Mdotn,Gamma,P0);

    % Obtaining average B and H values on the magnet boundary

    % {B} = B.n + {Btan}
    Bn = reconstruct(TdAM+TdAJ,Gamma,curl(NED)).*normals;
    % nx curlAxn from the outside
    Btano = cross(normals,reconstruct(TnAJ+TnAM,Gamma,DIV0),2);
    Btani = Btano + Mxn;
    avgB = Bn + 0.5*(Btano + Btani);
    
    % Computing the integral of (Mxn)x{B}
    
    for fieldID = 1:Nfields
        a = abc_alpha(fieldID,1);
        b = abc_alpha(fieldID,2);
        c = abc_alpha(fieldID,3);
        alpha = abc_alpha(fieldID,4);
        [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);
        Vels = Vel(X);
        shape_derivatives_mst(i,fieldID) = sum(W.* dot(Vels,cross(Mxn,avgB,2),2),1);
    end
%     forces_mst(i,:) = sum(W.* cross(Mxn,avgB,2),1)

%     [Vel,DVel] = getPolyVelDVel(1,1,1,1);

    shape_derivatives_bem(i,:) = PermanentMagnetShapeDerivativeBIEVP_dualnorm(Gamma,TnAJ_RWG+TnAM_RWG,TdAJ+TdAM,J,omega_src,mu,mu,M,abc_alpha);
    save("PMVP_Cuboid_5_dualnorm.mat","shape_derivatives_bem","shape_derivatives_mst","hvals");
end
