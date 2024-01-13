function [] = PMCFSP_forcesNtorques(meshfunction,vals)
    delete(gcp('nocreate'));
    funcInfo = functions(meshfunction);
    disp(['PMCFSP forces and torques invoked with: ', funcInfo.function]);
    disp("===================================================");
    format long;
    mu0 = 1;
    Nvals = size(vals,2);
    
    forces_mst = zeros(Nvals,3);
    forces_bem = forces_mst;
    torques_mst = forces_mst;
    torques_bem = forces_mst;
    hvals = 0*vals;

    rng(32);
    ra1 = rand(1,3);
    ra2 = rand(1,3);

    for i = 1:Nvals
        N = 2^vals(i);
        disp(N);
        %% SOLUTION DOMAIN
        bndmesh_i = meshfunction(N);

        % Bounding box
        bndmesh_e = mshSphere(N,9);
        bndmesh_e = bndmesh_e.translate([2 2 2]);
        assert(vecnorm(bndmesh_i.vtx-[2 2 2],2,2)<9);
        % Mesh size
        hvals(i) = sqrt(mean(bndmesh_i.ndv,1));
        
        Gamma_i = dom(bndmesh_i,3);
        normals_i = Gamma_i.qudNrm;
        
        %% Solving the transmission problem
        % Values that seem to work
        % M = @(X) ones(size(X,1),1)  * [10 5 1];
        % H0 = [3 2 1];
        M = @(X) ones(size(X,1),1)  * [10 0 0];
        H0 = [1 0 0];
        [psi_i,g_i,psi_e] = solveTpPMCFSP(bndmesh_i,bndmesh_e,mu0,M);

        %% Computing the force using equivalent charge model
        P1_i = fem(bndmesh_i,'P1');
        P0_i = fem(bndmesh_i,'P0');
        H0extended = repmat(H0,size(normals_i,1),1);

        [X_i,W_i] = Gamma_i.qud;

        Mvals = M(X_i);
        Mdotn = dot(Mvals,normals_i,2);

        psi_i_vals = reconstruct(psi_i,Gamma_i,P0_i);
        gradg_i_vals = reconstruct(g_i,Gamma_i,P1_i.grad);

        avgHtot = (-psi_i_vals-0.5 * Mdotn).*normals_i + gradg_i_vals + H0extended;

        fdensity = avgHtot .* Mdotn;
    
        forces_mst(i,:) = mu0 * sum(W_i.*fdensity,1)
    
        Xcg = [4 0 0];
        r = X_i-Xcg;
        torques_mst(i,:) = mu0 * sum(W_i.*cross(r,fdensity,2),1)
        
        %% Computing the BEM based force and torque
        
        [Vel1,DVel1] = getTransVelDVel([1 0 0]);
        [Vel2,DVel2] = getTransVelDVel([0 1 0]);
        [Vel3,DVel3] = getTransVelDVel([0 0 1]);
        
        f1 = SdBEMPMCFSP_ConstVEL(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel1,DVel1,mu0);
        f2 = SdBEMPMCFSP_ConstVEL(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel2,DVel2,mu0);
        f3 = SdBEMPMCFSP_ConstVEL(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel3,DVel3,mu0);

        forces_bem(i,:) = [f1 f2 f3]

        [Velr1,DVelr1] = getRotVelDVel([1 0 0],Xcg);
        [Velr2,DVelr2] = getRotVelDVel([0 1 0],Xcg);
        [Velr3,DVelr3] = getRotVelDVel([0 0 1],Xcg);
        
        t1 = SdBEMPMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Velr1,DVelr1,mu0,H0,M);
        t2 = SdBEMPMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Velr2,DVelr2,mu0,H0,M);
        t3 = SdBEMPMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Velr3,DVelr3,mu0,H0,M);

        torques_bem(i,:) = [t1 t2 t3]

        fname = "PMCFSP_forcesNtorques_" + funcInfo.function + ".mat";
        save(fname,"forces_mst","torques_mst","forces_bem","torques_bem","hvals");
    end

end