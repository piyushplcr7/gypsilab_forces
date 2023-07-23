% Panel oriented assembly for double integral with a kernel in 3D

% I: index subset of elements
% J: index subset of elements
% size(I, 1) must be equal to size(J, 1).

% Pairs of elements I x J share a vertex, edge or are identical.


function M = panel_assembly_VECTORIZED(mesh,kernel,trial_space,test_space,Intmat)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% FAR AWAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [I0,J0] = find(Intmat == 0);
    
    elts_I0 = mesh.elt(I0,:);
    elts_J0 = mesh.elt(J0,:);
    
    AI0 = mesh.vtx(elts_I0(:,1),:);
    BI0 = mesh.vtx(elts_I0(:,2),:);
    CI0 = mesh.vtx(elts_I0(:,3),:);
    
    AJ0 = mesh.vtx(elts_J0(:,1),:);
    BJ0 = mesh.vtx(elts_J0(:,2),:);
    CJ0 = mesh.vtx(elts_J0(:,3),:);
    
    ABCI0 = elts_I0;
    ABCJ0 = elts_J0;
    
    permI0 = repmat([1 2 3],size(elts_I0,1),1);
    permJ0 = permI0;

    [DCVxI0,DCVyI0] = GramianVectorized(AI0,BI0,CI0);
    [DCVxJ0,DCVyJ0] = GramianVectorized(AJ0,BJ0,CJ0);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% COMMON VERTEX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [I1,J1] = find(Intmat == 1);
    
    % Getting elements in matrix form corresponding to I1 and J1
    elts_I1 = mesh.elt(I1,:);
    elts_J1 = mesh.elt(J1,:);
    
    [intersection1,diffI1,diffJ1] = rowWiseIntersectionDiff(elts_I1,elts_J1);
    
    % Getting the vertices
    
    % Common vertex! AI1 = AJ1
    AI1 = mesh.vtx(intersection1,:);
    AJ1 = AI1;
    
    BI1 = mesh.vtx(diffI1(:,1),:);
    CI1 = mesh.vtx(diffI1(:,2),:);
    
    BJ1 = mesh.vtx(diffJ1(:,1),:);
    CJ1 = mesh.vtx(diffJ1(:,2),:);
    
    ABCI1 = [intersection1 diffI1];
    ABCJ1 = [intersection1 diffJ1];
    
    % Finding permutation, perm is such that elt(perm) = ABC
    permI1 = findPermVectorized(ABCI1,elts_I1);
    permJ1 = findPermVectorized(ABCJ1,elts_J1);

    [DCVxI1,DCVyI1] = GramianVectorized(AI1,BI1,CI1);
    [DCVxJ1,DCVyJ1] = GramianVectorized(AJ1,BJ1,CJ1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% COMMON EDGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [I2,J2] = find(Intmat == 2);
    
    elts_I2 = mesh.elt(I2,:);
    elts_J2 = mesh.elt(J2,:);
    
    [intersection2,diffI2,diffJ2] = rowWiseIntersectionDiff(elts_I2,elts_J2);
    
    % Getting the vertices
    
    % Common edge! AI2 = AJ2, BI2 = BJ2
    AI2 = mesh.vtx(intersection2(:,1),:);
    BI2 = mesh.vtx(intersection2(:,2),:);
    AJ2 = AI2;
    BJ2 = BI2;
    
    CI2 = mesh.vtx(diffI2,:);
    CJ2 = mesh.vtx(diffJ2,:);
    
    ABCI2 = [intersection2 diffI2];
    ABCJ2 = [intersection2 diffJ2];
    
    % Finding permutation, perm is such that elt(perm) = ABC
    permI2 = findPermVectorized(ABCI2,elts_I2);
    permJ2 = findPermVectorized(ABCJ2,elts_J2);

    [DCVxI2,DCVyI2] = GramianVectorized(AI2,BI2,CI2);
    [DCVxJ2,DCVyJ2] = GramianVectorized(AJ2,BJ2,CJ2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% IDENTICAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [I3,J3] = find(Intmat == 3);
    
    % Element I and J are identical
    elts_IJ3 = mesh.elt(I3,:);
    
    % Getting the vertices
    
    % Identical element! AI3 = AJ3, BI3 = BJ3, CI3 = CJ3
    AI3 = mesh.vtx(elts_IJ3(:,1),:);
    BI3 = mesh.vtx(elts_IJ3(:,2),:);
    CI3 = mesh.vtx(elts_IJ3(:,3),:);
    
    AJ3 = AI3;
    BJ3 = BI3;
    CJ3 = CI3;
    
    ABCI3 = elts_IJ3;
    ABCJ3 = elts_IJ3;
    
    % Finding permutation, perm is such that elt(perm) = ABC
    %permI3 = findPermVectorized(ABCI3,elts_IJ3);
    %permJ3 = findPermVectorized(ABCJ3,elts_IJ3);
    permI3 = repmat([1 2 3],size(elts_IJ3,1),1);
    permJ3 = permI3;

    [DCVxI3,DCVyI3] = GramianVectorized(AI3,BI3,CI3);
    DCVxJ3 = DCVxI3;
    DCVyJ3 = DCVyI3;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X, W] = quad4D(5); 
    Xss{1} = X{1}; Wss{1} = W{1}; % Identical
    Xss{2} = X{2}; Wss{2} = W{2}; % Common edge
    Xss{3} = X{3}; Wss{3} = W{3}; % Common vertex
    Xss{4} = X{4}; Wss{4} = W{4}; % far away

    % Making the matrix of quadrature points and weights for vectorization
    % 1st component of X points
    % X1Mat = [repmat(X{4}(:,1)',size(I0,1));... % Far away points
    %         repmat(X{3}(:,1)',size(I1,1));... % Common vertex points
    %         repmat(X{2}(:,1)',size(I2,1));... % Common edge points
    %         repmat(X{1}(:,1)',size(I3,1))]    % Identical points
    % 
    % X2Mat = [repmat(X{4}(:,2)',size(I0,1));... % Far away points
    %         repmat(X{3}(:,2)',size(I1,1));... % Common vertex points
    %         repmat(X{2}(:,2)',size(I2,1));... % Common edge points
    %         repmat(X{1}(:,2)',size(I3,1))]    % Identical points
    % 
    % Y1Mat = [repmat(X{4}(:,3)',size(I0,1));... % Far away points
    %         repmat(X{3}(:,3)',size(I1,1));... % Common vertex points
    %         repmat(X{2}(:,3)',size(I2,1));... % Common edge points
    %         repmat(X{1}(:,3)',size(I3,1))]    % Identical points
    % 
    % Y2Mat = [repmat(X{4}(:,4)',size(I0,1));... % Far away points
    %         repmat(X{3}(:,4)',size(I1,1));... % Common vertex points
    %         repmat(X{2}(:,4)',size(I2,1));... % Common edge points
    %         repmat(X{1}(:,4)',size(I3,1))]    % Identical points
    % 
    % WMat = [repmat(X{4}',size(I0,1));... % Far away points
    %         repmat(X{3}',size(I1,1));... % Common vertex points
    %         repmat(X{2}',size(I2,1));... % Common edge points
    %         repmat(X{1}',size(I3,1))]    % Identical points

    % Converting to the right reference element
    % % X4 = Xss{4}(:, 1:2);
    % % Y4 = Xss{4}(:, 3:4);
    % % X4(:,1) = X4(:,1) - X4(:,2);
    % % Y4(:,1) = Y4(:,1) - Y4(:,2);
    % X1Mat = X1Mat - X2Mat;
    % Y1Mat = Y1Mat - Y2Mat;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Computing all interactions 0 %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % No. of quadrature points for far away
    Nqud0 = size(X{4},1);
    % No. of quadrature points for common vertex
    Nqud1 = size(X{3},1);
    % No. of quadrature points for edge
    Nqud2 = size(X{2},1);
    % No. of quadrature points for identical elements
    Nqud3 = size(X{1},1);

    % X quadrature points
    XMat = [repmat([X{4}(:,1)-X{4}(:,2) X{4}(:,2)],size(I0,1),1);... % Nqud0 * NI0 X 2
            repmat([X{3}(:,1)-X{3}(:,2) X{3}(:,2)],size(I1,1),1);... % Nqud1 * NI1 X 2
            repmat([X{2}(:,1)-X{2}(:,2) X{2}(:,2)],size(I2,1),1);... % Nqud2 * NI2 X 2
            repmat([X{1}(:,1)-X{1}(:,2) X{1}(:,2)],size(I3,1),1)];   % Nqud3 * NI3 X 2

    YMat = [repmat([X{4}(:,3)-X{4}(:,4) X{4}(:,4)],size(I0,1),1);... % Nqud0 * NI0 X 2
            repmat([X{3}(:,3)-X{3}(:,4) X{3}(:,4)],size(I1,1),1);... % Nqud1 * NI1 X 2
            repmat([X{2}(:,3)-X{2}(:,4) X{2}(:,4)],size(I2,1),1);... % Nqud2 * NI2 X 2
            repmat([X{1}(:,3)-X{1}(:,4) X{1}(:,4)],size(I3,1),1)];   % Nqud3 * NI3 X 2

    WMat = [repmat(W{4},size(I0,1),1);... % Nqud0 *NI0 x 1
            repmat(W{3},size(I1,1),1);... % Nqud1 *NI1 x 1
            repmat(W{2},size(I2,1),1);... % Nqud2 *NI2 x 1
            repmat(W{1},size(I3,1),1)];   % Nqud3 *NI3 x 1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% chi_tau <-> i and xhat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    AIvec = [repelem(AI0,Nqud0,1);...
             repelem(AI1,Nqud1,1);...
             repelem(AI2,Nqud2,1);...
             repelem(AI3,Nqud3,1)];

    BIvec = [repelem(BI0,Nqud0,1);...
             repelem(BI1,Nqud1,1);...
             repelem(BI2,Nqud2,1);...
             repelem(BI3,Nqud3,1)];

    CIvec = [repelem(CI0,Nqud0,1);...
             repelem(CI1,Nqud1,1);...
             repelem(CI2,Nqud2,1);...
             repelem(CI3,Nqud3,1)];

    chi_tau = AIvec + (BIvec-AIvec).*XMat(:,1) + (CIvec-AIvec).*XMat(:,2);

%     AI0vec = repelem(AI0,Nqud,1);
%     BI0vec = repelem(BI0,Nqud,1);
%     CI0vec = repelem(CI0,Nqud,1);
%     chi_tau0 = AI0vec + (BI0vec-AI0vec).*XMat(:,1) + (CI0vec-AI0vec).*Xmat(:,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% chi_t <-> j and yhat %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    AJvec = [repelem(AJ0,Nqud0,1);...
             repelem(AJ1,Nqud1,1);...
             repelem(AJ2,Nqud2,1);...
             repelem(AJ3,Nqud3,1)];

    BJvec = [repelem(BJ0,Nqud0,1);...
             repelem(BJ1,Nqud1,1);...
             repelem(BJ2,Nqud2,1);...
             repelem(BJ3,Nqud3,1)];

    CJvec = [repelem(CJ0,Nqud0,1);...
             repelem(CJ1,Nqud1,1);...
             repelem(CJ2,Nqud2,1);...
             repelem(CJ3,Nqud3,1)];

    chi_t = AJvec + (BJvec-AJvec).*YMat(:,1) + (CJvec-AJvec).*YMat(:,2);

%     AJ0vec = repelem(AJ0,Nqud,1);
%     BJ0vec = repelem(BJ0,Nqud,1);
%     CJ0vec = repelem(CJ0,Nqud,1);
%     chi_t0 = AJ0vec + (BJ0vec-AJ0vec).*YMat(:,1) + (CJ0vec-AJ0vec).*YMat(:,2);

    % All permutations
    ABCI = [ABCI0; ABCI1; ABCI2; ABCI3];
    ABCJ = [ABCJ0; ABCJ1; ABCJ2; ABCJ3];

    % All permutations
    permI = [permI0; permI1; permI2; permI3];
    permJ = [permJ0; permJ1; permJ2; permJ3];

    Ivec = [I0;I1;I2;I3];
    Jvec = [J0;J1;J2;J3];
    NIvec = [size(I0,1);size(I1,1);size(I2,1);size(I3,1)];

    Nqudvec = [Nqud0;Nqud1;Nqud2;Nqud3];

    % All elements
    elts_I = mesh.elt(Ivec,:);
    elts_J = mesh.elt(Jvec,:);

    % DCVx, DCVy
    DCVxI = [repelem(DCVxI0,Nqud0,1);...
             repelem(DCVxI1,Nqud1,1);...
             repelem(DCVxI2,Nqud2,1);...
             repelem(DCVxI3,Nqud3,1)];

    DCVyI = [repelem(DCVyI0,Nqud0,1);...
             repelem(DCVyI1,Nqud1,1);...
             repelem(DCVyI2,Nqud2,1);...
             repelem(DCVyI3,Nqud3,1)];

    DCVxJ = [repelem(DCVxJ0,Nqud0,1);...
             repelem(DCVxJ1,Nqud1,1);...
             repelem(DCVxJ2,Nqud2,1);...
             repelem(DCVxJ3,Nqud3,1)];

    DCVyJ = [repelem(DCVyJ0,Nqud0,1);...
             repelem(DCVyJ1,Nqud1,1);...
             repelem(DCVyJ2,Nqud2,1);...
             repelem(DCVyJ3,Nqud3,1)];

    % Normals and area info from the mesh
    normals = mesh.nrm;
    vols = mesh.ndv;

    % Element info
    elts_I0 = mesh.elt(I0,:);
    elts_J0 = mesh.elt(J0,:);

    % Normals and area info at interacting pairs I0 and J0
%     vols_I0 = vols(I0);
%     vols_J0 = vols(J0);
%     g_tau_I0 = 2 * vols_I0; % I <-> x
%     g_t_J0 = 2 * vols_J0; % J <-> y
    
    % normals at I0 and J0
%     nrmx = normals(I0,:);
%     nrmy = normals(J0,:); 

    % Vectorizing for quadrature
%     g_tau_I0_vec = repelem(g_tau_I0,Nqud,1);
    g_tau_I_vec = 2 * [repelem(vols(I0),Nqud0,1);...
                       repelem(vols(I1),Nqud1,1);...
                       repelem(vols(I2),Nqud2,1);...
                       repelem(vols(I3),Nqud3,1)];

%     g_t_J0_vec = repelem(g_t_J0,Nqud,1);
    g_t_J_vec = 2 * [repelem(vols(J0),Nqud0,1);...
                     repelem(vols(J1),Nqud1,1);...
                     repelem(vols(J2),Nqud2,1);...
                     repelem(vols(J3),Nqud3,1)];

%     nrmx_vec = repelem(nrmx,Nqud,1);
    nrmx_vec = [repelem(normals(I0,:),Nqud0,1);...
                repelem(normals(I1,:),Nqud1,1);...
                repelem(normals(I2,:),Nqud2,1);...
                repelem(normals(I3,:),Nqud3,1)];

%     nrmy_vec = repelem(nrmy,Nqud,1);
    nrmy_vec = [repelem(normals(J0,:),Nqud0,1);...
                repelem(normals(J1,:),Nqud1,1);...
                repelem(normals(J2,:),Nqud2,1);...
                repelem(normals(J3,:),Nqud3,1)];

    % Vectorized kernel evaluation 
%     KerVec0 = kernel(chi_tau0,chi_t0,chi_t0-chi_tau0);
    KerVec = kernel(chi_tau,chi_t,chi_t-chi_tau);

    
    % Info about trial and test spaces
    [~,elt2dof_tr] = trial_space.dof;
    [~,elt2dof_ts] = test_space.dof;
    Qtr = size(elt2dof_tr,2);
    Qts = size(elt2dof_ts,2);
    
    rsf_tr = trial_space.rsf;
    rsf_ts = test_space.rsf;
    tr_typ = trial_space.typ;
    ts_typ = test_space.typ;

    tr_opr = trial_space.opr;
    ts_opr = test_space.opr;

    % Output matrix
    M = zeros(test_space.ndof,trial_space.ndof);

    switch tr_typ
        case 'P0' % P0 trial function and P0 test, 1 RSF each
            Psix = rsf_ts{1}(1) .* g_tau_I_vec;
            Psiy = rsf_tr{1}(1) .* g_t_J_vec;
            
            % Quadrature before summation, contains NI0,1,2,3 interactions
            % at the same time
            integrands = WMat.*Psix.*KerVec.*Psiy; 
            integrals = reduceIntegrand(integrands,NIvec,Nqudvec);
            
            %local_matrix(ii,jj) = dot(Wh,Psix .* Ker .* Psiy);
%             M(Ivec, Jvec) = M(Ivec, Jvec) + integrals;
            % The above operation doesn't work as expected, need to use
            % linear indices
            indices = sub2ind(size(M),Ivec,Jvec);
            M(indices) = M(indices) + integrals;
        
        case 'P1' % P1 trial function
            switch tr_opr
                case '[psi]' % Trial: P1, Test: P1, 3 RSFs each
                    Psix = g_tau_I_vec.*[rsf_ts{1}(XMat) rsf_ts{2}(XMat) rsf_ts{3}(XMat)];
                    Psiy = g_t_J_vec.*[rsf_tr{1}(YMat) rsf_tr{2}(YMat) rsf_tr{3}(YMat)];

                    % Injecting weights and kernel into Psix
                    Psix = Psix.*WMat.*KerVec;
                    
                    % If Psix = [a1 a2 a3], 
                    % Psix_triple_vector = [a1 a2 a3
                    %                       a1 a2 a3
                    %                       a1 a2 a3]
                    Psix_triple_vector = repmat(Psix,3,1);

                    % If Psiy = [b1 b2 b3]
                    % Psiy_triple_row = [b1 b2 b3 b1 b2 b3 b1 b2 b3];
                    Psiy_triple_row = repmat(Psiy,1,3);

                    % Reshaping to get long vectors
                    % Psix_triple_vector = [a1 
                    %                       a1
                    %                       a1
                    %                       a2
                    %                       a2
                    %                       a2
                    %                       a3
                    %                       a3
                    %                       a3]
                    Psix_triple_vector = reshape(Psix_triple_vector,9*size(XMat,1),1);

                    % Psiy_triple_row    = [b1 
                    %                       b2
                    %                       b3
                    %                       b1
                    %                       b2
                    %                       b3
                    %                       b1
                    %                       b2
                    %                       b3]
                    Psiy_triple_row = reshape(Psiy_triple_row,9*size(YMat,1),1);

                    % All i,j combinations
                    % Combinations =       [a1 b1
                    %                       a1 b2
                    %                       a1 b3
                    %                       a2 b1
                    %                       a2 b2
                    %                       a2 b3
                    %                       a3 b1 
                    %                       a3 b2
                    %                       a3 b3]
                    Psix_triple_vector = reshape(Psix_triple_vector,9*size(XMat,1),1);
                    combinations = Psix_triple_vector .* Psiy_triple_row;

                    % Reducing the combinations
                    % Combinations =
                    % [a1b1 a1b2 a1b3 a2b1 a2b2 a2b3 a3b1 a3b2 a3b3]
                    combinations = reshape(combinations,size(XMat,1),9);

                    % Summing up for quadrature
                    combinations = reduceIntegrand(combinations,NIvec,Nqudvec);

                    % Filling the global matrix
                    % Getting the right indices
%                     M( reshape(repmat(ABCI,3,1),9*size(ABCI,1),1) ,...
%                        repmat(reshape(ABCJ,3*size(ABCJ,1)),3,1) ) = ...
%                     M( reshape(repmat(ABCI,3,1),9*size(ABCI,1),1) ,...
%                        repmat(reshape(ABCJ,3*size(ABCJ,1)),3,1) ) +...
%                     combinations(:);
                    indices = sub2ind(size(M),...
                            reshape(repmat(ABCI,3,1),9*size(ABCI,1),1),...
                            repmat(reshape(ABCJ,3*size(ABCJ,1),1),3,1));
                    M(indices) = M(indices) + combinations(:);
                    
                case 'grad[psi]' % Trial gradP1, test: nxRWG, both have 3 RSF
                          assert(strcmp(ts_typ,'RWG'));
                          assert(strcmp(ts_opr,'nx[psi]'));

                          dofs_I = elt2dof_ts(Ivec,:);

                          for ii = 1:Qts
                                ip1 = mod(permI(:,ii),3)+1;
                                ip2 = mod(ip1,3)+1;

                                % Fluxes
                                flux0 = (2*(mesh.elt(I0,ip1) < mesh.elt(I0,ip2))-1);
                                flux1 = (2*(mesh.elt(I1,ip1) < mesh.elt(I1,ip2))-1);
                                flux2 = (2*(mesh.elt(I2,ip1) < mesh.elt(I2,ip2))-1);
                                flux3 = (2*(mesh.elt(I3,ip1) < mesh.elt(I3,ip2))-1);

                                flux = [repelem(flux0,Nqud0,1);...
                                        repelem(flux1,Nqud1,1);...
                                        repelem(flux2,Nqud2,1);...
                                        repelem(flux3,Nqud3,1)];
                                
                                
                                % RTO reference element
                                Psix_ref = flux.*[rsf_ts{ii}{1}(XMat)...
                                                                  rsf_ts{ii}{2}(XMat)];

                                % Transforming the basis, g_tau/(2 volsi)
                                % cancels
                                Psix = (BIvec-AIvec).*Psix_ref(:,1)...
                                    + (CIvec-AIvec).*Psix_ref(:,2);
                                
                                % evaluating nxPsi
                                nxPsix = cross(nrmx_vec,Psix,2);

                                for jj = 1:Qtr
                                    % Reference gradient
                                    dPsiy_ref = [rsf_tr{jj}{1} rsf_tr{jj}{2}];

                                    % Transformed element obtained by 
                                    % multiplication with DCV
                                    dPsiy = DCVxJ * dPsiy_ref(1)...
                                           +DCVyJ * dPsiy_ref(2);

                                    % adding g_t
                                    dPsiy = dPsiy.*g_t_J_vec;
                                    
                                    integrand_ii_jj = WMat.*KerVec.*dot(nxPsix,dPsiy,2);
                                    integral_ii_jj = reduceIntegrand(integrand_ii_jj,NIvec,Nqudvec);
                                    
                                    dofsI_at_permI = evalRowWisePerm(dofs_I,permI);
                                    M(dofsI_at_permI(:,ii),ABCJ(:,jj)) = M(dofsI_at_permI(:,ii),ABCJ(:,jj)) + integral_ii_jj;

                                end

                          end

                case 'nxgrad[psi]'  % Trial,test: nxgrad(P1), 3 RSFs
                    % Panel pairs we want to check
                    itest = 1;
                    jtest = 5;
                    % Find this pair in Ivec, Jvec
                    fi = find(Ivec == itest);
                    fj = find(Jvec(fi) == jtest);
                    found_idx = fi(fj);
                    assert(Ivec(found_idx) == itest);
                    assert(Jvec(found_idx) == jtest);

                    % Checking DCV
                    found_indices = findIndicesInBigVec(itest,jtest,I0,I1,I2,I3,J0,J1,J2,J3,Nqudvec);
                    found_indices = found_indices(1)

%                     Ai=AIvec(found_indices,:)
%                     Bi=BIvec(found_indices,:)
%                     Ci=CIvec(found_indices,:)

%                     Aj=AJvec(found_indices,:)
%                     Bj=BJvec(found_indices,:)
%                     Cj=CJvec(found_indices,:)
                    
%                     DCVj = [DCVxJ(found_indices,:)' DCVyJ(found_indices,:)']
%                     DCVi = [DCVxI(found_indices,:)' DCVyI(found_indices,:)']
                    

                    for ii = 1:Qts
                        % Reference gradient
                        dPsix_ref = [rsf_ts{ii}{1} rsf_ts{ii}{2}];

                        % Transformed element obtained by 
                        % multiplication with DCV
                        dPsix = DCVxI * dPsix_ref(1)...
                               +DCVyI * dPsix_ref(2);

                        % adding g_tau
                        dPsix = dPsix.*g_tau_I_vec;

                        nxdPsix = cross(nrmx_vec,dPsix,2);
                       
                        for jj = 1:Qtr
                            % Reference gradient
                            dPsiy_ref = [rsf_tr{jj}{1} rsf_tr{jj}{2}];

                            % Transformed element obtained by 
                            % multiplication with DCV
                            dPsiy = DCVxJ * dPsiy_ref(1)...
                                   +DCVyJ * dPsiy_ref(2);

                            % adding g_t
                            dPsiy = dPsiy.*g_t_J_vec;

                            nxdPsiy = cross(nrmy_vec,dPsiy,2);

                            integrand_ii_jj = WMat.*KerVec.*dot(nxdPsix,nxdPsiy,2);
                            integral_ii_jj = reduceIntegrand(integrand_ii_jj,NIvec,Nqudvec);
                            indices = sub2ind(size(M),ABCI(:,ii),ABCJ(:,jj));
%                             M(ABCI(:,ii),ABCJ(:,jj)) = M(ABCI(:,ii),ABCJ(:,jj)) + integral_ii_jj;
%                             M(indices) = M(indices) + integral_ii_jj;

                            for var = 1:size(ABCI,1)
                                M(ABCI(var,ii),ABCJ(var,jj)) = M(ABCI(var,ii),ABCJ(var,jj)) + integral_ii_jj(var); 

                            end
%                             Ntif = size(testidxfind,1);
%                             if Ntif ~= 0
%                                 for yo = 1:Ntif
%                                  fprintf("Interaction of panels %d and %d for RSFs %d and %d gives %.9f \n", Ivec(found_idx),Jvec(found_idx),ii,jj,integral_ii_jj(found_idx));
%                                 end
%                             end

                            % Galerkin matrix entries we want to check
%                             it = 1; jt = 1;
%                             ijidx = sub2ind(size(M),it,jt);
%                             findijidx = find(indices = ididx);

                        end

                    end

                case 'n*[psi]'
                    switch ts_opr

                        case '[psi]' % Trial: ntimes(P1), Test: P0
                            % Kernel has to be vectorial! ker(x,y).n(y)
                            Kerdotn_vec = dot(KerVec,nrmy_vec,2);

                            for ii = 1:Qts
                                Psix = rsf_ts{ii}(XMat).*g_tau_I_vec;
                                for jj = 1:Qtr
                                    Psiy = rsf_tr{jj}(YMat).*g_t_J_vec;

                                    integrand_ii_jj = WMat.*Kerdotn_vec.*Psix.*Psiy;
                                    integral_ii_jj = reduceIntegrand(integrand_ii_jj,NIvec,Nqudvec);
%                                     M(Ivec,ABCJ(:,jj)) = M(Ivec,ABCJ(:,jj)) + integral_ii_jj;
                                    for var = 1:size(Ivec,1)
                                        M(Ivec(var),ABCJ(var,jj)) = M(Ivec(var),ABCJ(var,jj)) + integral_ii_jj(var);
                                    end
                                end
                            end

                        case 'n*[psi]' % Trial: ntimes(P1), Test: ntimes(P1)
                            ndotn = dot(nrmx_vec,nrmy_vec);

                            for ii = 1:Qts
                                Psix = rsf_ts{ii}(XMat).*g_tau_I_vec;
                                for jj = 1:Qtr
                                    Psiy = rsf_tr{jj}(YMat).*g_t_J_vec;
                                    integrand_ii_jj = WMat.*KerVec.*Psix.*Psiy.*ndotn;
                                    integral_ii_jj = reduceIntegrand(integrand_ii_jj,NIvec,Nqudvec);
                                    M(ABCI(:,ii),ABCJ(:,jj)) = M(ABCI(:,ii),ABCJ(:,jj)) + integral_ii_jj;
                                end
                            end
 


                    end

            end % Switch tr_opr

        case 'RWG'
            switch tr_opr
                case '[psi]'

                    dofs_I = elt2dof_ts(Ivec,:);
                    dofs_J = elt2dof_tr(Jvec,:);

                    for ii = 1:Qts
                        ip1 = mod(permI(:,ii),3)+1;
                        ip2 = mod(ip1,3)+1;

                        % Fluxes
                        flux0i = (2*(mesh.elt(I0,ip1) < mesh.elt(I0,ip2))-1);
                        flux1i = (2*(mesh.elt(I1,ip1) < mesh.elt(I1,ip2))-1);
                        flux2i = (2*(mesh.elt(I2,ip1) < mesh.elt(I2,ip2))-1);
                        flux3i = (2*(mesh.elt(I3,ip1) < mesh.elt(I3,ip2))-1);

                        fluxi = [repelem(flux0i,Nqud0,1);...
                                repelem(flux1i,Nqud1,1);...
                                repelem(flux2i,Nqud2,1);...
                                repelem(flux3i,Nqud3,1)];
                        
                        
                        % RTO reference element
                        Psix_ref = fluxi.*[rsf_ts{ii}{1}(XMat)...
                                                          rsf_ts{ii}{2}(XMat)];

                        % Transforming the basis, g_tau/(2 volsi)
                        % cancels
                        Psix = (BIvec-AIvec).*Psix_ref(:,1)...
                            + (CIvec-AIvec).*Psix_ref(:,2);

                        for jj = 1:Qtr
                            jp1 = mod(permJ(:,jj),3)+1;
                            jp2 = mod(jp1,3)+1;
    
                            % Fluxes
                            flux0j = (2*(mesh.elt(J0,jp1) < mesh.elt(J0,jp2))-1);
                            flux1j = (2*(mesh.elt(J1,jp1) < mesh.elt(J1,jp2))-1);
                            flux2j = (2*(mesh.elt(J2,jp1) < mesh.elt(J2,jp2))-1);
                            flux3j = (2*(mesh.elt(J3,jp1) < mesh.elt(J3,jp2))-1);
    
                            fluxj = [repelem(flux0j,Nqud0,1);...
                                    repelem(flux1j,Nqud1,1);...
                                    repelem(flux2j,Nqud2,1);...
                                    repelem(flux3j,Nqud3,1)];
                            
                            
                            % RTO reference element
                            Psiy_ref = fluxj.*[rsf_tr{jj}{1}(YMat)...
                                                              rsf_tr{jj}{2}(YMat)];
    
                            % Transforming the basis, g_tau/(2 volsi)
                            % cancels
                            Psiy = (BJvec-AJvec).*Psiy_ref(:,1)...
                                + (CJvec-AJvec).*Psiy_ref(:,2);

                            if size(KerVec,2) == 1
                                integrand_ii_jj = WMat.*KerVec.*dot(Psix,Psiy,2);
                            elseif size(KerVec,2) == 3
                                % Need to perform {kernel x trial}.test
                                kerxtrial_vec = cross(KerVec,Psiy,2);
                                integrand_ii_jj = WMat.*dot(kerxtrial_vec,Psix,2);
                            end

                            integral_ii_jj = reduceIntegrand(integrand_ii_jj,NIvec,Nqudvec);
                            
                            dofsI_at_permI = evalRowWisePerm(dofs_I,permI);
                            dofsJ_at_permJ = evalRowWisePerm(dofs_J,permJ);
                            M(dofsI_at_permI(:,ii),dofsJ_at_permJ(:,jj)) = M(dofsI_at_permI(:,ii),dofsJ_at_permJ(:,jj)) + integral_ii_jj;
                        end
                    end

            end


        
    end % Switch tr_typ


end % Function end
   