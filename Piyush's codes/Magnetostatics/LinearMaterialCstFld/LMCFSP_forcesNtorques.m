function [] = LMCFSP_forcesNtorques(meshfunction,vals)
    
    funcInfo = functions(meshfunction);
    disp(['LMCFSP forces and torques invoked with: ', funcInfo.function]);
    disp("===================================================");
    format long;
    mu = 4;
    mu0 = 2;
    Nvals = size(vals,2);
    
    forces_mst = zeros(Nvals,3);
    forces_bem = forces_mst;
    torques_mst = forces_mst;
    torques_bem = forces_mst;
    hvals = 0*vals;

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
        H0 = [10 3 1];
        [psi_i,g_i,psi_e] = solveTPLMCFSP(bndmesh_i,bndmesh_e,mu,mu0,H0);

        %% Computing the force using MST formula
        % Reconstructing Bn and Ht
        P1_i = fem(bndmesh_i,'P1');
        P0_i = fem(bndmesh_i,'P0');
        H0extended = repmat(H0,size(normals_i,1),1);
        H0t = cross(normals_i,cross(H0extended,normals_i,2),2);
        Ht = reconstruct(g_i,Gamma_i,P1_i.grad) + H0t;
    
        Bn = -mu0 * reconstruct(psi_i,Gamma_i,P0_i) + mu0 * dot(H0extended,normals_i,2);
    
        forces_mst(i,:) = ForceMstTP(Gamma_i,Bn,vecnorm(Ht,2,2),mu0,mu);
    
        Xcg = [4 0 0];
        torques_mst(i,:) = TorqueMstTP(Gamma_i,Bn,vecnorm(Ht,2,2),mu0,mu,Xcg);
        
        %% Computing the BEM based force and torque
        
        [Vel1,DVel1] = getTransVelDVel([1 0 0]);
        [Vel2,DVel2] = getTransVelDVel([0 1 0]);
        [Vel3,DVel3] = getTransVelDVel([0 0 1]);
    
        f1 = SdBEMLMCFSP_ConstVEL(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel1,DVel1,mu0,mu,H0);
        f2 = SdBEMLMCFSP_ConstVEL(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel2,DVel2,mu0,mu,H0);
        f3 = SdBEMLMCFSP_ConstVEL(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Vel3,DVel3,mu0,mu,H0);

        forces_bem(i,:) = [f1 f2 f3];

        [Velr1,DVelr1] = getRotVelDVel([1 0 0],Xcg);
        [Velr2,DVelr2] = getRotVelDVel([0 1 0],Xcg);
        [Velr3,DVelr3] = getRotVelDVel([0 0 1],Xcg);

        t1 = SdBEMLMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Velr1,DVelr1,mu0,mu,H0);
        t2 = SdBEMLMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Velr2,DVelr2,mu0,mu,H0);
        t3 = SdBEMLMCFSP(bndmesh_i,bndmesh_e,psi_i,g_i,psi_e,Velr3,DVelr3,mu0,mu,H0);
    
        torques_bem(i,:) = [t1 t2 t3];

        fname = "LMCFSP_forcesNtorques_" + funcInfo.function + ".mat";
        save(fname,"forces_mst","torques_mst","forces_bem","torques_bem","hvals");
    end

end