function [] = LMCFSP_ALT_dualnormCuboid(vals)
    delete(gcp('nocreate'));
    disp('LMCFSP ALT dualnorm invoked with Cuboid ');
    disp("===================================================");

    format long;
    mu = 4;
    mu0 = 2;
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
        bndmesh_i = getMeshCuboid5Translated(floor(N/3),[1 0.5 1]);

        % Bounding box
        bndmesh_e = mshSphere(40*floor(N^0.7),4);
        
        % Mesh size
        hvals(i) = sqrt(mean(bndmesh_i.ndv,1));

        order = 3;
        Gamma_i = dom(bndmesh_i,order);
        normals_i = Gamma_i.qudNrm;
        
        %% Solving the transmission problem
        H0 = [10 3 1];

        % Full solution (Dirichlet excitation)
        [psi_i,g_i,psi_e] = solveTPLMCFSP_ALT(bndmesh_i,bndmesh_e,mu,mu0,H0,order);
        
        % Constructing traces for the full solution
        P1_i = fem(bndmesh_i,'P1');
        P0_i = fem(bndmesh_i,'P0');
        [X_i,W_i] = Gamma_i.qud;

        % Reconstructing Bn and Ht
        gradudotn = -reconstruct(psi_i,Gamma_i,P0_i);
        sgradu = reconstruct(g_i,Gamma_i,grad(P1_i));
    
        Ht = vecnorm(sgradu,2,2);
        Bn = mu0 * gradudotn;
        
        jump_mu_inv = 1/mu0 - 1/mu;
        jump_mu = mu0 - mu;
        
        fdensity = 0.5 * ((Bn).^2*jump_mu_inv - (Ht).^2*jump_mu).* normals_i;
        
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
        shape_derivatives_bem(i,:) = SdBEMLMCFSP_dualnorm_ALT(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,mu0,mu,H0,abc_alpha,order);
%         disp("BEM shape derivatives computed!");

        fname = "LMCFSP_dualnorm_ALT_getMeshCuboid5.mat";
        save(fname,"shape_derivatives_bem","shape_derivatives_mst","hvals");
    end

end