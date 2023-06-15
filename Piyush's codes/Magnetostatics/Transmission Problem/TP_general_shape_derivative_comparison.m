% New transmission problem script

addpath(genpath("../../../"));
clear; clc; close all;
format long;
% (mui+mue)/(mui-mue)
mu = 4;
mu0 = 2;
vals = 5:12;
Nvals = size(vals,2);
sd_vol = zeros(Nvals,1);
sd_bem = sd_vol;
hvals = 0*vals;

for i = 1:Nvals
    N = 2^vals(i);
    disp(N);
    %% SOLUTION DOMAIN
    % Cube size and position
    L = 2*[1 1 1];
    T = [5 5 3];
    
    % Cube domain
    %bndmesh = bndmeshCubeTranslated(N,L,T);
    
    % Spherical domain
    bndmesh = mshSphere(N,1);
    bndmesh = bndmesh.translate(T);
    
%     mesh = mshCube(N,L);
%     mesh = mesh.translate(T);
%     %mesh = mesh.sub(1);
%     bndmesh = mesh.bnd;
    
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
    
    %% Computing Volume based Shape Derivative
    
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

    a = str2num(getenv("AA"));
    b = str2num(getenv("BB"));
    c = str2num(getenv("CC"));
    k = str2num(getenv("KK"));
    
    [Vel,DVel] = getPolyVelDVel(a,b,c,k);

    sd_vol(i) = ShapeDervTpVol(Gamma,Bn,Ht,mu0,mu,Vel)

    %% Computing BEM shape derivative

    % Projecting Psi to RWG
    Psi_RWG = proj(Psivals,Gamma,RWG);

    sd_bem(i) = SdBemTPVP(bndmesh,Psi_RWG,g,J,omega_src,Vel,DVel,mu0,mu)

    fname = "sd_sph_poly_";
    suffix = ".mat";
    fname = strcat(fname,int2str(a),int2str(b),int2str(c),int2str(k),suffix);

    save(fname,"sd_vol","sd_bem");
end
