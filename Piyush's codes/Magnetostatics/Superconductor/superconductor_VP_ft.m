function [] = superconductor_VP_ft(meshfunction,vals)
    funcInfo = functions(meshfunction);
    disp(['SC VP forces and torques invoked with: ', funcInfo.function]);
    disp("===================================================");
    format long;
    
    Nvals = size(vals,2);
    hvals = vals;

    forces_mst = zeros(Nvals,3);
    torques_mst = forces_mst;
    forces_bem = forces_mst;
    torques_bem = torques_mst;

    for i = 1:Nvals
        N = 2^vals(i);
        disp(N)
        %for N=50
        %% SOLUTION DOMAIN
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
        
        %% Solving the problem and obtaining the Neumann trace
        [TnA0,TnA] = solve_superconductor(bndmesh,J,mesh_src);
        
        %% Plotting the computed B field
        
        %plot_field(TnA,bndmesh,J,omega_src);
        
        %% Computing forces
        % Coefficients for zero Dirichlet Trace
        TdA = TnA * 0;
        forces_mst(i,:) = MstForceFromA(TdA,TnA,Gamma)'
        
        % Shape Derivative computation of force
        % Translation fields
        Nux = @(X) ones(size(X,1),1)*[1 0 0];
        Nuy = @(X) ones(size(X,1),1)*[0 1 0];
        Nuz = @(X) ones(size(X,1),1)*[0 0 1];
        
        sd_e1 = -SuperConductorShapeDerivativeT3(bndmesh,TnA,Nux,omega_src,J);
        sd_e2 = -SuperConductorShapeDerivativeT3(bndmesh,TnA,Nuy,omega_src,J);
        sd_e3 = -SuperConductorShapeDerivativeT3(bndmesh,TnA,Nuz,omega_src,J);
    
        forces_bem(i,:) = [sd_e1 sd_e2 sd_e3]
        
        %% Computing torques
        Xcg = [4 0 0];
        torques_mst(i,:) = MstTorqueFromA(TdA,TnA,Gamma,Xcg)'
    
        % Shape Derivative computation of torque
        % Getting Rotational Vels and DVels
        [Velxr,DVelxr] = getRotVelDVel([1 0 0],Xcg);
        [Velyr,DVelyr] = getRotVelDVel([0 1 0],Xcg);
        [Velzr,DVelzr] = getRotVelDVel([0 0 1],Xcg);
    
    
        sdt_e1 = SuperConductorShapeDerivative(bndmesh,TnA,Velxr,DVelxr,omega_src,J)
        sdt_e2 = SuperConductorShapeDerivative(bndmesh,TnA,Velyr,DVelyr,omega_src,J)
        sdt_e3 = SuperConductorShapeDerivative(bndmesh,TnA,Velzr,DVelzr,omega_src,J)
    
        torques_bem(i,:) = [sdt_e1 sdt_e2 sdt_e3]
         % torques_bem(i,:) = ptorque'
    
        fname = "SCVP_ft_" + funcInfo.function + ".mat";
        save(fname,"forces_bem","forces_mst", "torques_bem","torques_mst","hvals");
        % save("SC_VP_Sph.mat","forces_mst","forces_bem","torques_mst","torques_bem","hvals");

    end

end