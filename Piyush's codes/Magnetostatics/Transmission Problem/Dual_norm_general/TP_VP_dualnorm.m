function [] = TP_VP_dualnorm(meshfunction,vals)
    funcInfo = functions(meshfunction);
    disp(['TP VP dualnorm invoked with: ', funcInfo.function]);
    disp("===================================================");

    % New transmission problem script
    
%     clear; clc; close all;
    format long;

    mu = 4;
    mu0 = 2;
    %vals = 5:12;
    Nvals = size(vals,2);

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
    
    % Number of fields
    Nfields = size(abc_alpha,1);
  
    hvals = 0*vals;
    hvals_src = hvals;
    
    for i = 1:Nvals
        N = 2^vals(i);
        disp(N);
        %% SOLUTION DOMAIN
        % Cube size and position
        bndmesh = meshfunction(N);
        
        % Mesh size
        hvals(i) = sqrt(mean(bndmesh.ndv,1));
        
        Gamma = dom(bndmesh,3);
        normals = Gamma.qudNrm;
        
        %% Source
        if strcmp(funcInfo.function,'getMeshCube')
            disp('getMeshCube Nsrc');
            N_src = floor(8*N^0.7);
        elseif strcmp(funcInfo.function,'getMeshSphere')
            disp('getMeshSphere Nsrc');
            N_src = floor(3*N);
        elseif strcmp(funcInfo.function,'getMeshTetraNew')
            disp('getMeshTetraNew Nsrc');
            N_src = floor((N^2)/350);
        end

        R0 = 2;
        r0 = .5;
        [J,mesh_src] = get_torus_source(N_src,R0,r0);
        omega_src = dom(mesh_src,3);
        hvals_src(i) = sqrt(mean(mesh_src.ndv,1));
        
        %% Solving the transmission problem
        % These are traces from the exterior
        % Psi lies in the space nxgradP1 and g lies in the space NED
        [Psi,g] = solveTransmissionProblem(bndmesh,J,omega_src,mu,mu0);
%         disp("Solution computed!");

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
        
%         disp("MST shape derivatives computed!");
    
        %% Computing forces and torques using BEM shape derivative
    
        % Projecting Psi to RWG
        Psi_RWG = proj(Psivals,Gamma,RWG);
        shape_derivatives_bem(i,:) = SdBemTPVP_dualnorm(bndmesh,Psi_RWG,g,J,omega_src,mu0,mu,abc_alpha);
%         disp("BEM shape derivatives computed!");

        fname = "TPVP_dualnorm_" + funcInfo.function + ".mat";
        save(fname,"shape_derivatives_bem","shape_derivatives_mst","hvals");
    end

end
