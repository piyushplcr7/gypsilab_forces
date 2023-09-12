function [] = superconductor_SP_dualnorm(meshfunction,vals)

    funcInfo = functions(meshfunction);
    disp(['SC SP dualnorm invoked with: ', funcInfo.function]);
    disp("===================================================");

    format long;
    mui = 1;
    mue = 1;
%     vals = 5:12;
    Nvals = size(vals,2);
    hvals = 0*vals;
    
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
        %N = vals(i);
    
        disp(N)
        %for N=50
        %% SOLUTION DOMAIN
        % Cube size and position
        bndmesh = meshfunction(N);
        
        % Mesh size
        hvals(i) = sqrt(mean(bndmesh.ndv,1));
    
        % Dom object
        Gamma = dom(bndmesh,3);
        
        %% Current source
        N_src = N;
        R0 = 2;
        r0 = .5;
        [J,mesh_src] = get_torus_source(N_src,R0,r0);
        omega_src = dom(mesh_src,3);
        
        %% Solving the problem and obtaining the Dirichlet Trace of scalar pot
        [Tdu,Tnu] = solve_superconductor_scalar_potential(bndmesh,J,mesh_src);
        
        %% Plotting the computed B field
        
    %     plot(bndmesh);
    %     hold on;
    %     plot_field_scalar_potential(bndmesh,Tdu,Tnu,J,omega_src);
        %plot_field(TnA,bndmesh,J,omega_src);
    
        %% Computing the volume based shape derivative
    
        for fieldID = 1:Nfields
            a = abc_alpha(fieldID,1);
            b = abc_alpha(fieldID,2);
            c = abc_alpha(fieldID,3);
            alpha = abc_alpha(fieldID,4);
            [Vel,DVel] = getCosVelDVel(a,b,c,alpha+1);
            shape_derivatives_mst(i,fieldID) = ScSd_SP_Vol(bndmesh,Tdu,Tnu,J,omega_src,Vel,DVel);
            
        end
    
        shape_derivatives_bem(i,:) = ScSd_SP_BEM_dualnorm(bndmesh,Tdu,Tnu,J,omega_src,abc_alpha);
        
        fname = "SCSP_dualnorm_" + funcInfo.function + ".mat";
        save(fname,"shape_derivatives_bem","shape_derivatives_mst","hvals");
    end

end