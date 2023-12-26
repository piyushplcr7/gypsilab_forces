function [] = superconductor_SP_ft(meshfunction,vals)
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
        %N = vals(i);
    
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
        
        %% Solving the problem and obtaining the Dirichlet Trace of scalar pot
        [Tdu,Tnu] = solve_superconductor_scalar_potential(bndmesh,J,mesh_src);
        
        %% Plotting the computed B field
        
    %     plot(bndmesh);
    %     hold on;
    %     plot_field_scalar_potential(bndmesh,Tdu,Tnu,J,omega_src);
        %plot_field(TnA,bndmesh,J,omega_src);
    
        %% Computing the volume based shape derivative
        % Computing forces
        [Vel1,DVel1] = getTransVelDVel([1 0 0]);
        [Vel2,DVel2] = getTransVelDVel([0 1 0]);
        [Vel3,DVel3] = getTransVelDVel([0 0 1]);
    
        f1 = ScSd_SP_Vol(bndmesh,Tdu,Tnu,J,omega_src,Vel1,DVel1);
        f2 = ScSd_SP_Vol(bndmesh,Tdu,Tnu,J,omega_src,Vel2,DVel2);
        f3 = ScSd_SP_Vol(bndmesh,Tdu,Tnu,J,omega_src,Vel3,DVel3);
    
        forces_mst(i,:) = [f1 f2 f3]
    
        % Computing torques
        Xcg = [4 0 0];
        [Velr1,DVelr1] = getRotVelDVel([1 0 0],Xcg);
        [Velr2,DVelr2] = getRotVelDVel([0 1 0],Xcg);
        [Velr3,DVelr3] = getRotVelDVel([0 0 1],Xcg);
    
        t1 = ScSd_SP_Vol(bndmesh,Tdu,Tnu,J,omega_src,Velr1,DVelr1);
        t2 = ScSd_SP_Vol(bndmesh,Tdu,Tnu,J,omega_src,Velr2,DVelr2);
        t3 = ScSd_SP_Vol(bndmesh,Tdu,Tnu,J,omega_src,Velr3,DVelr3);
    
        torques_mst(i,:) = [t1 t2 t3]
    
        %% Computing BEM based shape derivative
        % Computing forces
        fbem1 = ScSdF_SP_BEM(bndmesh,Tdu,Tnu,J,omega_src,Vel1);
        fbem2 = ScSdF_SP_BEM(bndmesh,Tdu,Tnu,J,omega_src,Vel2);
        fbem3 = ScSdF_SP_BEM(bndmesh,Tdu,Tnu,J,omega_src,Vel3);
    
        forces_bem(i,:) = [fbem1 fbem2 fbem3]
    
        %Computing torques
        tbem1 = ScSd_SP_BEM(bndmesh,Tdu,Tnu,J,omega_src,Velr1,DVelr1)
        tbem2 = ScSd_SP_BEM(bndmesh,Tdu,Tnu,J,omega_src,Velr2,DVelr2)
        tbem3 = ScSd_SP_BEM(bndmesh,Tdu,Tnu,J,omega_src,Velr3,DVelr3)
        
        torques_bem(i,:) = [tbem1 tbem2 tbem3]
    
        fname = "SCSP_ft_" + funcInfo.function + ".mat";
        save(fname,"forces_bem","forces_mst", "torques_bem","torques_mst","hvals");
        % save("SC_SP_Sph.mat","forces_mst","torques_mst","forces_bem","torques_bem","hvals");
    end


end