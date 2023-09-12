function [] = PM_SP_dualnorm_gradphi(meshfunction,vals,M)
    funcInfo = functions(meshfunction);
    MfuncInfo = functions(M);
    disp(['PM SP dualnorm invoked with: ', funcInfo.function, ' ', MfuncInfo.function]);
    disp("===================================================");
    
    % Script for constant Magnetization
    format long;
    
    mu = 1;
    mu0 = mu;
%     vals = 5:13;
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
    
    for i = 1:Nvals
        N = 2^vals(i);
        disp(N);
        %% SOLUTION DOMAIN
        % Cube size and position
        bndmesh = meshfunction(N);
        
        % Mesh size
        hvals(i) = sqrt(mean(bndmesh.ndv,1));
        
        Gamma = dom(bndmesh,12);
        normals = Gamma.qudNrm;
    
        % Function to determine the interior of the magnet
        %interior = @(X) (vecnorm(X-T,2,2) < 1);
        
        %% Source
        N_src = N;
        R0 = 2;
        r0 = .5;
        [J_orig,mesh_src] = get_torus_source(N_src,R0,r0);
        omega_src = dom(mesh_src,3);
        
        %% Solving the Magnet problem (SP)
        % Constant Magnetization function
%         M = @(X) ones(size(X,1),1) * [1 0 0];
        % Modifying J source
        J = @(x,y,z) 1 * J_orig(x,y,z);
    
        % Solving for exterior traces
        [psi,g] = solveMagnetProblemSimplifiedSP(Gamma,mu,mu,J,omega_src,M);
    
        %% Computing the B field for verification
    %     plot_field_magnet(TdAM,TnAM_RWG,bndmesh,J,omega_src,mu0,interior);
    %     figure;
    %     plot_field_magnet(TdAJ,TnAJ_RWG,bndmesh,J,omega_src,mu0,interior);
    
        %% Computing forces and torques from Equivalent charge model
        
        [X,W] = Gamma.qud;
        Mvals = M(X);
        Mdotn = dot(Mvals,normals,2);
    
        % Obtaining H values on the magnet boundary
        % {H} = {H.n} + Htan
        HJ = compute_vecpot_curl(J,omega_src,X);
    
        % Htan = grad_T g
        P1 = fem(bndmesh,'P1');
        Htan = reconstruct(g,Gamma,P1.grad);
    
        % Normal component is stored in the Neumann trace
        P0 = fem(bndmesh,'P0');
        Hnout = reconstruct(psi,Gamma,P0);
        Hnin = Hnout-Mdotn;
        avgHn = 0.5*(Hnout+Hnin);
        avgH = avgHn.*normals + Htan + HJ;
        
        % Computing the surface integral of (M.n) {H}
    %     forces_mst(i,:) = mu0 * sum(W.* Mdotn .* avgH,1)
    
        for fieldID = 1:Nfields
            a = abc_alpha(fieldID,1);
            b = abc_alpha(fieldID,2);
            c = abc_alpha(fieldID,3);
            alpha = abc_alpha(fieldID,4);
            [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);
            Vels = Vel(X);
            shape_derivatives_mst(i,fieldID) = sum(W.* Mdotn .* dot(avgH,Vels,2),1);
        end
    
        shape_derivatives_bem(i,:) = PermanentMagnetShapeDerivativeBIESP_dualnorm(Gamma,g,psi,J,omega_src,mu0,M,abc_alpha);

        fname = "PMSP_dualnorm_" + funcInfo.function + "_" + MfuncInfo.function + ".mat";
        save(fname,"shape_derivatives_bem","shape_derivatives_mst","hvals");
    end

end