function [] = DielSC_dualnorm(meshfunction,vals)
    
    funcInfo = functions(meshfunction);
    disp(['DielSC dualnorm invoked with: ', funcInfo.function]);
    disp("===================================================");

    % New transmission problem script
%     clear; clc; close all;
    format long;
    % (mui+mue)/(mui-mue)
    epsilon = 4;
    epsilon0 = 2;
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
        mesh_src = mshSphere(N,1);
        % mesh_src = mshCube(N,L);
        mesh_src = mesh_src.translate([5 4 3]);
        omega_src = dom(mesh_src,3);
        % Constant surface charge density
        rho = @(X) ones(size(X,1),1);
        
        %% Solving the transmission problem 
        [Tnu,Tdu] = solveTpDielSourceCharge(bndmesh,epsilon,epsilon0,rho,omega_src);
        
        %% Computing the MST based force and torque

        Gamma = dom(bndmesh,3);
        P0 = fem(bndmesh,'P0');
        P1 = fem(bndmesh,'P1');
    
        jumpepsilon = epsilon0 - epsilon;
        jumpInvepsilon = 1/epsilon0 - 1/epsilon;
    
        grad_g_vals = reconstruct(Tdu,Gamma,grad(P1));
        alpha_vals = epsilon0 * reconstruct(Tnu,Gamma,P0);
        [X,W] = Gamma.qud;
        normals = Gamma.qudNrm;
    
        fdensity = 0.5 * (jumpepsilon * vecnorm(grad_g_vals,2,2).^2 - jumpInvepsilon * alpha_vals.^2 ) .* normals;
        
        for fieldID = 1:Nfields
            a = abc_alpha(fieldID,1);
            b = abc_alpha(fieldID,2);
            c = abc_alpha(fieldID,3);
            alpha = abc_alpha(fieldID,4);
            [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);

            Vels = Vel(X);
            shape_derivatives_mst(i,fieldID) = sum(W.*dot(fdensity,Vels,2),1);
        end

%         disp("MST shape derivatives computed!");
    
        % BEM based shape derivative
        shape_derivatives_bem(i,:) = sdBEMTpDielSC_dualnorm(bndmesh,epsilon,epsilon0,Tnu,Tdu,rho,omega_src,abc_alpha);
%         disp("BEM shape derivatives computed!");

        fname = "DielSC_dualnorm_" + funcInfo.function + ".mat";
        save(fname,"shape_derivatives_bem","shape_derivatives_mst","hvals");
    end

end