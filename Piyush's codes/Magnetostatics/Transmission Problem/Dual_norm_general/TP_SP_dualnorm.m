function [] = TP_SP_dualnorm(meshfunction,vals)
    
    funcInfo = functions(meshfunction);
    disp(['TP SP dualnorm invoked with: ', funcInfo.function]);
    disp("===================================================");

    % New transmission problem script
%     clear; clc; close all;
    format long;
    % (mui+mue)/(mui-mue)
    mu = 4;
    mu0 = 2;
    %vals = 5:12;
    Nvals = size(vals,2);
    hvals = vals;
    
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
        N_src = N;
        R0 = 2;
        r0 = .5;
        [J,mesh_src] = get_torus_source(N_src,R0,r0);
        omega_src = dom(mesh_src,3);
        
        %% Solving the transmission problem using scalar potential formulation
        % These are traces from the exterior
        % Psi lies in the space nxgradP1 and g lies in the space NED
        [Tnu,Tdu] = solveTpScalPotBIE(bndmesh,mu,mu0,J,omega_src);
%         disp("Solution computed!");
        
        %% Computing the MST based force and torque
    
        for fieldID = 1:Nfields
            a = abc_alpha(fieldID,1);
            b = abc_alpha(fieldID,2);
            c = abc_alpha(fieldID,3);
            alpha = abc_alpha(fieldID,4);
            [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);
            shape_derivatives_mst(i,fieldID) = shapeDer_ScalPot_Vol_TP(bndmesh,mu,mu0,Tnu,Tdu,J,omega_src,Vel,DVel);
        end

%         disp("MST shape derivatives computed!");
    
        % BEM based shape derivative
        shape_derivatives_bem(i,:) = shapDervTranPrbScalPotBIE_dualnorm(bndmesh,Tdu,Tnu,J,omega_src,mu0,mu,abc_alpha);
%         disp("BEM shape derivatives computed!");

        fname = "TPSP_dualnorm_" + funcInfo.function + ".mat";
        save(fname,"shape_derivatives_bem","shape_derivatives_mst","hvals");
    end

end