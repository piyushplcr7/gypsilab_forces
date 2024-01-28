function [] = LMCFVP_dualnorm(meshfunction,vals)
    
    funcInfo = functions(meshfunction);
    disp(['LMCFVP dualnorm invoked with: ', funcInfo.function]);
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
        bndmesh_i = meshfunction(N);

        % Bounding box
        bndmesh_e = mshSphere(N,9);
        bndmesh_e = bndmesh_e.translate([2 2 2]);
        
        % Mesh size
        hvals(i) = sqrt(mean(bndmesh_i.ndv,1));
        
        Gamma_i = dom(bndmesh_i,3);
        Gamma_e = dom(bndmesh_e,3);
        normals_i = Gamma_i.qudNrm;
        
        %% Solving the transmission problem
        % These are traces from the exterior
        B_0 = [1 0 0];
        [Psi_i,g_i,Psi_e] = solveTPLMCFVP(bndmesh_i,bndmesh_e,mu,mu0,B_0);

        % Reconstructing Bn and Ht
        NED_i = fem(bndmesh_i,'NED'); 
        P1_i = fem(bndmesh_i,'P1');
        DIV0_i = nxgrad(P1_i); 
        RWG_i = fem(bndmesh_i,'RWG');
    
        % Bn = curlA.n = curlTg
        Bn = reconstruct(g_i,Gamma_i,NED_i.curl);
        Btotn = normals_i * B_0' + Bn;
    
        Psivals_i = reconstruct(Psi_i,Gamma_i,DIV0_i);
        B0_tan = B_0 - (normals_i * B_0').*normals_i;
        % -ve sign for psivals_i because of chosen convention in derivation
        Htot_tan = 1/mu0 * (B0_tan - cross(normals_i,Psivals_i,2));
        Htot_tan = vecnorm(Htot_tan,2,2);

        jump_mu_inv = 1/mu0 - 1/mu;
        jump_mu = mu0 - mu;
        [X_i,W_i] = Gamma_i.qud;
        fdensity = 0.5 * ((Btotn).^2*jump_mu_inv - (Htot_tan).^2*jump_mu).* normals_i;

        % Projecting traces to RWG Spaces
        RWG_e = fem(bndmesh_e,'RWG');
        P1_e = fem(bndmesh_e,'P1');
        DIV0_e = nxgrad(P1_e);
        Psivals_e = reconstruct(Psi_e,Gamma_e,DIV0_e);
        Psie_RWG = proj(Psivals_e,Gamma_e,RWG_e);
        Psii_RWG = proj(Psivals_i,Gamma_i,RWG_i);
        
        %% Computing the MST based force and torque
        
        for fieldID = 1:Nfields
            a = abc_alpha(fieldID,1);
            b = abc_alpha(fieldID,2);
            c = abc_alpha(fieldID,3);
            alpha = abc_alpha(fieldID,4);
            [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);

            Vels = Vel(X_i);
            shape_derivatives_mst(i,fieldID) = sum(W_i.*dot(fdensity,Vels,2),1);
        end

%         disp("MST shape derivatives computed!");
        
        % BEM based shape derivative
        shape_derivatives_bem(i,:) = SdBemLMCFVP_dualnorm(bndmesh_i,bndmesh_e,Psii_RWG,g_i,Psie_RWG,mu0,mu,B_0,abc_alpha);
%         disp("BEM shape derivatives computed!");

        fname = "LMCFVP_dualnorm_" + funcInfo.function + ".mat";
        save(fname,"shape_derivatives_bem","shape_derivatives_mst","hvals");
    end

end