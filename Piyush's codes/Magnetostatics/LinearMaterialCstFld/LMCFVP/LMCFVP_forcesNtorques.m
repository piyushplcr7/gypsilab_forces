function [] = LMCFVP_forcesNtorques(meshfunction,vals)
    delete(gcp('nocreate'));
    funcInfo = functions(meshfunction);
    disp(['LMCFVP forces and torques invoked with: ', funcInfo.function]);
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
        [bndmesh_i,bndmesh_e] = meshfunction(N);

        % Bounding box
        % bndmesh_e = mshSphere(N,9);
        % bndmesh_e = bndmesh_e.translate([2 2 2]);
        
        % Mesh size
        hvals(i) = sqrt(mean(bndmesh_i.ndv,1));
        
        Gamma_i = dom(bndmesh_i,3);
        Gamma_e = dom(bndmesh_e,3);
        normals_i = Gamma_i.qudNrm;
        
        %% Solving the transmission problem
        % These are traces from the exterior
        B_0 = [1 0 0];
        [Psi_i,g_i,Psi_e] = solveTPLMCFVP(bndmesh_i,bndmesh_e,mu,mu0,B_0);

        %% MST Based force and torque on bndmesh_i
        % Force computation
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
    
        forces_mst(i,:) = ForceMstTP(Gamma_i,Btotn,Htot_tan,mu0,mu)
    
        % Torque computation
        Xcg = [4 0 0];
        torques_mst(i,:) = TorqueMstTP(Gamma_i,Btotn,Htot_tan,mu0,mu,Xcg)
        
        %% BEM Based SD Computation
        % Projecting traces to RWG Spaces
        RWG_e = fem(bndmesh_e,'RWG');
        P1_e = fem(bndmesh_e,'P1');
        DIV0_e = nxgrad(P1_e);
        Psivals_e = reconstruct(Psi_e,Gamma_e,DIV0_e);
        Psie_RWG = proj(Psivals_e,Gamma_e,RWG_e);
        Psii_RWG = proj(Psivals_i,Gamma_i,RWG_i);
    
        [Vel1,DVel1] = getTransVelDVel([1 0 0]);
        [Vel2,DVel2] = getTransVelDVel([0 1 0]);
        [Vel3,DVel3] = getTransVelDVel([0 0 1]);
        
        ff1 = SdBemLMCFVP_ConstVel(bndmesh_i,bndmesh_e,Psii_RWG,g_i,Psie_RWG,Vel1,DVel1,mu0,mu,B_0);
        ff2 = SdBemLMCFVP_ConstVel(bndmesh_i,bndmesh_e,Psii_RWG,g_i,Psie_RWG,Vel2,DVel2,mu0,mu,B_0);
        ff3 = SdBemLMCFVP_ConstVel(bndmesh_i,bndmesh_e,Psii_RWG,g_i,Psie_RWG,Vel3,DVel3,mu0,mu,B_0);
    
        forces_bem(i,:) = [ff1 ff2 ff3]
    
        [Velr1,DVelr1] = getRotVelDVel([1 0 0],Xcg);
        [Velr2,DVelr2] = getRotVelDVel([0 1 0],Xcg);
        [Velr3,DVelr3] = getRotVelDVel([0 0 1],Xcg);
    
        t1 = SdBemLMCFVP(bndmesh_i,bndmesh_e,Psii_RWG,g_i,Psie_RWG,Velr1,DVelr1,mu0,mu,B_0);
        t2 = SdBemLMCFVP(bndmesh_i,bndmesh_e,Psii_RWG,g_i,Psie_RWG,Velr2,DVelr2,mu0,mu,B_0);
        t3 = SdBemLMCFVP(bndmesh_i,bndmesh_e,Psii_RWG,g_i,Psie_RWG,Velr3,DVelr3,mu0,mu,B_0);
    
        torques_bem(i,:) = [t1 t2 t3]

        fname = "LMCFVP_forcesNtorques_" + funcInfo.function + ".mat";
        save(fname,"forces_mst","torques_mst","forces_bem","torques_bem","hvals");
    end

end