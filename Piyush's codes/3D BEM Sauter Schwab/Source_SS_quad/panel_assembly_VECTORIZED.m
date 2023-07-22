% Panel oriented assembly for double integral with a kernel in 3D

% I: index subset of elements
% J: index subset of elements
% size(I, 1) must be equal to size(J, 1).

% Pairs of elements I x J share a vertex, edge or are identical.


function M = panel_assembly_VECTORIZED(mesh,kernel,trial_space,test_space, I, J,Intmat)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% FAR AWAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [I0,J0] = find(Intmat == 0);
    
    elts_I0 = bndmesh.elt(I0,:);
    elts_J0 = bndmesh.elt(J0,:);
    
    AI0 = bndmesh.vtx(elts_I0(1,:),:);
    BI0 = bndmesh.vtx(elts_I0(2,:),:);
    CI0 = bndmesh.vtx(elts_I0(3,:),:);
    
    AJ0 = bndmesh.vtx(elts_J0(1,:),:);
    BJ0 = bndmesh.vtx(elts_J0(2,:),:);
    CJ0 = bndmesh.vtx(elts_J0(3,:),:);
    
    ABCI0 = elts_I0;
    ABCJ0 = elts_J0;
    
    permI0 = repmat([1 2 3],size(elts_IJ0,1),1);
    permJ0 = permI0;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% COMMON VERTEX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [I1 J1] = find(Intmat == 1);
    
    % Getting elements in matrix form corresponding to I1 and J1
    elts_I1 = bndmesh.elt(I1,:);
    elts_J1 = bndmesh.elt(J1,:);
    
    % cell_elts_I1 = mat2cell(elts_I1,ones(size(I1,1),1));
    % cell_elts_J1 = mat2cell(elts_J1,ones(size(J1,1),1));
    % 
    % % Apply intersection element wise
    % intersection1 = cellfun(@intersect,cell_elts_I1,cell_elts_J1,'UniformOutput',false);
    % % Apply setdiff element wise
    % diffI1 = cellfun(@setdiff,cell_elts_I1,intersection1,'UniformOutput',false);
    % diffJ1 = cellfun(@setdiff,cell_elts_J1,intersection1,'UniformOutput',false);
    % 
    % % Converting back to matrix
    % intersection1 = cell2mat(intersection1);
    % diffI1 = cell2mat(diffI1);
    % diffI2 = cell2mat(diffI2);
    
    [intersection1,diffI1,diffJ1] = rowWiseIntersectionDiff(elts_I1,elts_J1);
    
    % Getting the vertices
    
    % Common vertex! AI1 = AJ1
    AI1 = bndmesh.vtx(intersection1,:);
    
    BI1 = bndmesh.vtx(diffI1(:,1),:);
    CI1 = bndmesh.vtx(diffI1(:,2),:);
    
    BJ1 = bndmesh.vtx(diffJ1(:,1),:);
    CJ1 = bndmesh.vtx(diffJ1(:,2),:);
    
    ABCI1 = [intersection1 diffI1];
    ABCJ1 = [intersection1 diffJ1];
    
    % Finding permutation, perm is such that elt(perm) = ABC
    permI1 = findPermVectorized(ABCI1,elts_I1);
    permJ1 = findPermVectorized(ABCJ1,elts_J1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% COMMON EDGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [I2,J2] = find(Intmat == 2);
    
    elts_I2 = bndmesh.elt(I2,:);
    elts_J2 = bndmesh.elt(J2,:);
    
    [intersection2,diffI2,diffJ2] = rowWiseIntersectionDiff(elts_I2,elts_J2);
    
    % Getting the vertices
    
    % Common edge! AI2 = AJ2, BI2 = BJ2
    AI2 = bndmesh.vtx(intersection2(1,:),:);
    BI2 = bndmesh.vtx(intersection2(2,:),:);
    
    CI2 = bndmesh.vtx(diffI2,:);
    CJ2 = bndmesh.vtx(diffJ2,:);
    
    ABCI2 = [intersection2 diffI2];
    ABCJ2 = [intersection2 diffJ2];
    
    % Finding permutation, perm is such that elt(perm) = ABC
    permI2 = findPermVectorized(ABCI2,elts_I2);
    permJ2 = findPermVectorized(ABCJ2,elts_J2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% IDENTICAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [I3,J3] = find(Intmat == 3);
    
    % Element I and J are identical
    elts_IJ3 = bndmesh.elt(I3,:);
    
    % Getting the vertices
    
    % Identical element! AI3 = AJ3, BI3 = BJ3, CI3 = CJ3
    AIJ3 = bndmesh.vtx(elts_IJ3(1,:),:);
    BIJ3 = bndmesh.vtx(elts_IJ3(2,:),:);
    CIJ3 = bndmesh.vtx(elts_IJ3(3,:),:);
    
    ABCI3 = elts_IJ3;
    ABCJ3 = elts_IJ3;
    
    % Finding permutation, perm is such that elt(perm) = ABC
    %permI3 = findPermVectorized(ABCI3,elts_IJ3);
    %permJ3 = findPermVectorized(ABCJ3,elts_IJ3);
    permI3 = repmat([1 2 3],size(elts_IJ3,1),1);
    permJ3 = permI3;
    
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


    % Computing all interactions 0
    Nqud = size(X{4},1);
    XMat = repmat([X{4}(:,1)-X{4}(:,2) X{4}(:,2)],size(I0,1),1);
    YMat = repmat([X{4}(:,3)-X{4}(:,4) X{4}(:,4)],size(I0,1),1);
    WMat = repmat(W{4},size(I0,1));

    % chi_tau <-> i and xhat
    % chi_t <-> j and yhat
    AI0vec = repelem(AI0,Nqud,1);
    BI0vec = repelem(BI0,Nqud,1);
    CI0vec = repelem(CIO,Nqud,1);
    chi_tau0 = AI0vec + (BI0vec-AI0vec).*XMat(:,1) + (CI0vec-AI0vec).*Xmat(:,2);

    AJ0vec = repelem(AJ0,Nqud,1);
    BJ0vec = repelem(BJ0,Nqud,1);
    CJ0vec = repelem(CJ0,Nqud,1);
    chi_t0 = AJ0vec + (BJ0vec-AJ0vec).*YMat(:,1) + (CJ0vec-AJ0vec).*YMat(:,2);

    KerVec0 = kernel(chi_tau0,chi_t0,chi_t0-chi_tau0);

    
    
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    M = zeros(test_space.ndof,trial_space.ndof);
    
     %M = sparse(test_space.ndof, trial_space.ndof);
    
    % Vector storing volume of the mesh elements
    vols = mesh.ndv;
    
    normals = mesh.nrm;
    

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
    
    % Convention: panel i <-> chi_tau //////<-> hbasisx
    %             panel j <-> chi_t   //////<-> hbasisy
    
    % Double loop over the mesh elements
    L = size(I, 1);
    
    for elt = 1:L
        i = I(elt);
        j = J(elt);
        
        
        dofs_i = elt2dof_ts(i, :);
        dofs_j = elt2dof_tr(j, :);
        
        area_i = vols(i);
        
        nrm = normals(j,:);
        
        nrmx = normals(i,:);
        % The Jacobian of transformation from the panel in 3D to reference
        % triangle of area 0.5
        g_tau = 2 * area_i;
        
%        for j = J'
          area_j = vols(j);
          % The Jacobian of transformation from the panel in 3D to reference
          % triangle of area 0.5
          g_t = 2 * area_j;
          
          %fprintf("Reached element i,j = %d, %d",i,j);
          
          % Finding the relation between panels i and j 
          intersection = intersect(mesh.elt(i,:),mesh.elt(j,:));
          vtcs = mesh.vtx(intersection,:);
          l = length(intersection);
          switch l
              case 0
                  relation = "far_away";
                  % Vertices for elt i
                  Ai = mesh.vtx(mesh.elt(i,1),:)';
                  Bi = mesh.vtx(mesh.elt(i,2),:)';
                  Ci = mesh.vtx(mesh.elt(i,3),:)';
                  % Vertices for elt j
                  Aj = mesh.vtx(mesh.elt(j,1),:)';
                  Bj = mesh.vtx(mesh.elt(j,2),:)';
                  Cj = mesh.vtx(mesh.elt(j,3),:)';
                  % Parameterizations
                  chi_tau = @(xhat) Ai + [Bi-Ai Ci-Ai]*xhat;
                  chi_t = @(yhat) Aj + [Bj-Aj Cj-Aj]*yhat;
                  ABC_elti = mesh.elt(i,:);
                  ABC_eltj = mesh.elt(j,:);
                  
                  
                  perm_i = 1:3;
                  perm_j = 1:3;
                  
                  % E is the Jacobian Matrix of the mapping for panel j
                  E = [Bj-Aj Cj-Aj]; 
                  % Gram matrix
                  EtE = E' * E; 
                  Dxy = inv(EtE);
                  
                  DCVx = Dxy(1, 1) * E(:, 1) + Dxy(1, 2) * E(:, 2);
                  DCVy = Dxy(2, 1) * E(:, 1) + Dxy(2, 2) * E(:, 2);
                  
                  DCV = [DCVx DCVy];

                  % Ei is the Jacobian Matrix of the Mapping for panel i
                  Ei = [Bi-Ai Ci-Ai]; 
                  % Gram Matrix
                  EitEi = Ei' * Ei; Dxy = inv(EitEi);
                  
                  DCVx = Dxy(1, 1) * Ei(:, 1) + Dxy(1, 2) * Ei(:, 2);
                  DCVy = Dxy(2, 1) * Ei(:, 1) + Dxy(2, 2) * Ei(:, 2);
                  
                  DCVi = [DCVx DCVy];
                  
                  % This is done to transform the reference element
                  % from (0,0), (1,0), (1,1) to
                  % (0,0), (1,0) , (0,1)
                  Xh = Xss{4}(:, 1:2);
                  Yh = Xss{4}(:, 3:4);
                  
                  Xh(:,1) = Xh(:,1) - Xh(:,2);
                  Yh(:,1) = Yh(:,1) - Yh(:,2);
                  Wh = Wss{4};
                  
                  
              case 1
                  relation = "common_vertex";
                  Ai = vtcs(1,:)'; Aj = Ai;
                  diffi = setdiff(mesh.elt(i,:),intersection);
                  BCi = mesh.vtx(diffi,:);
                  Bi = BCi(1,:)';
                  Ci = BCi(2,:)';
                  diffj = setdiff(mesh.elt(j,:),intersection);
                  BCj = mesh.vtx(diffj,:);
                  Bj = BCj(1,:)';
                  Cj = BCj(2,:)';
                  ABC_elti = [intersection diffi(1) diffi(2)];
                  ABC_eltj = [intersection diffj(1) diffj(2)];
                  % 0,0 for common point
                  chi_tau = @(xhat) Ai + [Bi-Ai Ci-Ai]*xhat;
                  chi_t = @(yhat) Aj + [Bj-Aj Cj-Aj]*yhat;
                  
                  
                  perm_i = [find(mesh.elt(i, :) == ABC_elti(1)), ...
                            find(mesh.elt(i, :) == ABC_elti(2)), ...
                            find(mesh.elt(i, :) == ABC_elti(3))];
                  
                  perm_j = [find(mesh.elt(j, :) == ABC_eltj(1)), ...
                            find(mesh.elt(j, :) == ABC_eltj(2)), ...
                            find(mesh.elt(j, :) == ABC_eltj(3))];
                  
                  
                  E = [Bj-Aj Cj-Aj]; EtE = E' * E; Dxy = inv(EtE);
                  
                  DCVx = Dxy(1, 1) * E(:, 1) + Dxy(1, 2) * E(:, 2);
                  DCVy = Dxy(2, 1) * E(:, 1) + Dxy(2, 2) * E(:, 2);
                  
                  DCV = [DCVx DCVy];
                  
                  Ei = [Bi-Ai Ci-Ai]; EitEi = Ei' * Ei; Dxy = inv(EitEi);
                  
                  DCVx = Dxy(1, 1) * Ei(:, 1) + Dxy(1, 2) * Ei(:, 2);
                  DCVy = Dxy(2, 1) * Ei(:, 1) + Dxy(2, 2) * Ei(:, 2);
                  
                  DCVi = [DCVx DCVy];
                  
                  
                  
                  Xh = Xss{3}(:, 1:2);
                  Yh = Xss{3}(:, 3:4);
                  
                  Xh(:,1) = Xh(:,1) - Xh(:,2);
                  Yh(:,1) = Yh(:,1) - Yh(:,2);
                  Wh = Wss{3};
                  
                  
              case 2
                  relation = "common_edge";
                  Ai = vtcs(1,:)'; Aj = Ai;
                  Bi = vtcs(2,:)'; Bj = Bi;
                  ci = setdiff(mesh.elt(i,:),intersection);
                  Ci = mesh.vtx(ci,:)';
                  cj = setdiff(mesh.elt(j,:),intersection);
                  Cj = mesh.vtx(cj,:)';
                  ABC_elti = [intersection(1) intersection(2) ci];
                  ABC_eltj = [intersection(1) intersection(2) cj];
                  % 0,0 and 1,0 for common points
                  chi_tau = @(xhat) Ai + [Bi-Ai Ci-Ai]*xhat;
                  chi_t = @(yhat) Aj + [Bj-Aj Cj-Aj]*yhat;
                  
                  perm_i = [find(mesh.elt(i, :) == ABC_elti(1)), ...
                            find(mesh.elt(i, :) == ABC_elti(2)), ...
                            find(mesh.elt(i, :) == ABC_elti(3))];
                  
                  perm_j = [find(mesh.elt(j, :) == ABC_eltj(1)), ...
                            find(mesh.elt(j, :) == ABC_eltj(2)), ...
                            find(mesh.elt(j, :) == ABC_eltj(3))];
                  
                  E = [Bj-Aj Cj-Aj]; EtE = E' * E; Dxy = inv(EtE);
                  
                  DCVx = Dxy(1, 1) * E(:, 1) + Dxy(1, 2) * E(:, 2);
                  DCVy = Dxy(2, 1) * E(:, 1) + Dxy(2, 2) * E(:, 2);
                  
                  DCV = [DCVx DCVy];
                  
                  Ei = [Bi-Ai Ci-Ai]; EitEi = Ei' * Ei; Dxy = inv(EitEi);
                  
                  DCVx = Dxy(1, 1) * Ei(:, 1) + Dxy(1, 2) * Ei(:, 2);
                  DCVy = Dxy(2, 1) * Ei(:, 1) + Dxy(2, 2) * Ei(:, 2);
                  
                  DCVi = [DCVx DCVy];
                  
                  
                  
                  Xh = Xss{2}(:, 1:2);
                  Yh = Xss{2}(:, 3:4);
                  
                  Xh(:,1) = Xh(:,1) - Xh(:,2);
                  Yh(:,1) = Yh(:,1) - Yh(:,2);
                  Wh = Wss{2};
                  
              case 3
                  relation = "identical";
                  Aj = vtcs(1,:)'; Ai = Aj;
                  Bj = vtcs(2,:)'; Bi = Bj;
                  Cj = vtcs(3,:)'; Ci = Cj;
                  chi_tau = @(xhat) Ai + [Bi-Ai Ci-Ai]*xhat;
                  chi_t = @(yhat) Aj + [Bj-Aj Cj-Aj]*yhat;
                  ABC_elti = intersection;
                  ABC_eltj = intersection;
                  
                  perm_i = [find(mesh.elt(i, :) == ABC_elti(1)), ...
                            find(mesh.elt(i, :) == ABC_elti(2)), ...
                            find(mesh.elt(i, :) == ABC_elti(3))];
                  
                  perm_j = [find(mesh.elt(j, :) == ABC_eltj(1)), ...
                            find(mesh.elt(j, :) == ABC_eltj(2)), ...
                            find(mesh.elt(j, :) == ABC_eltj(3))];
                  
                  E = [Bj-Aj Cj-Aj]; EtE = E' * E; Dxy = inv(EtE);
                  
                  DCVx = Dxy(1, 1) * E(:, 1) + Dxy(1, 2) * E(:, 2);
                  DCVy = Dxy(2, 1) * E(:, 1) + Dxy(2, 2) * E(:, 2);
                  
                  DCV = [DCVx DCVy];
                  
                  Ei = [Bi-Ai Ci-Ai]; EitEi = Ei' * Ei; Dxy = inv(EitEi);
                  
                  DCVx = Dxy(1, 1) * Ei(:, 1) + Dxy(1, 2) * Ei(:, 2);
                  DCVy = Dxy(2, 1) * Ei(:, 1) + Dxy(2, 2) * Ei(:, 2);
                  
                  DCVi = [DCVx DCVy];
                  
                  
                  
                  Xh = Xss{1}(:, 1:2);
                  Yh = Xss{1}(:, 3:4);
                  
                  
                  Xh(:,1) = Xh(:,1) - Xh(:,2);
                  Yh(:,1) = Yh(:,1) - Yh(:,2);
                  Wh = Wss{1};
          end
          
          
          local_matrix = zeros(Qts,Qtr);

          
          
          if iscell(kernel)
            Ker{1} = kernel{1}(chi_tau(Xh')',chi_t(Yh')',chi_t(Yh')'-chi_tau(Xh')');
            Ker{2} = kernel{2}(chi_tau(Xh')',chi_t(Yh')',chi_t(Yh')'-chi_tau(Xh')');
            Ker{3} = kernel{3}(chi_tau(Xh')',chi_t(Yh')',chi_t(Yh')'-chi_tau(Xh')');
          
          else
              
            Ker = kernel(chi_tau(Xh')',chi_t(Yh')',chi_t(Yh')'-chi_tau(Xh')');
          end
          
          switch tr_typ
              
              
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%% P0 FINITE ELEMENTS BASIS FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              case 'P0' 
          
                  for ii = 1:Qts
                      Psix = rsf_ts{ii}(Xh) .* g_tau;

                      for jj = 1:Qtr
                        Psiy = rsf_tr{jj}(Yh) .* g_t;
                        local_matrix(ii,jj) = dot(Wh,Psix .* Ker .* Psiy);
                      end

                  end
                  M(i, j) = M(i, j) + local_matrix;
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%%% P1 FINITE ELEMENTS BASIS FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              case 'P1'
                  
                  switch tr_opr

                      case '[psi]' % Trial: P1, Test: P1
                          
                      for ii = 1:Qts
                          Psix = rsf_ts{ii}(Xh) .* g_tau;

                          for jj = 1:Qtr
                            Psiy = rsf_tr{jj}(Yh) .* g_t;
                            local_matrix(ii,jj) = dot(Wh,Psix .* Ker .* Psiy);
                          end

                      end
                      M(ABC_elti, ABC_eltj) = M(ABC_elti, ABC_eltj) + local_matrix;

                      case 'grad[psi]' % Trial gradP1, test: nxRWG
                          assert(strcmp(ts_typ,'RWG'));
                          assert(strcmp(ts_opr,'nx[psi]'));

                          volsi = vols(i);
                          elti  = mesh.elt(i, :);

                          for ii = 1:Qts
                                ip1 = mod(perm_i(ii),3)+1;
                                ip2 = mod(ip1,3)+1;
                                
                                % Flux through the face
                                flux = (2*(elti(ip1) < elti(ip2))-1);
                                
                                dPsi1i = rsf_ts{ii}{1}(Xh); % Sizes N X 1
                                dPsi2i = rsf_ts{ii}{2}(Xh);
                                
                                % Size 2 X N
                                diP  = flux * [dPsi1i dPsi2i]'; % RT0 reference element 2X1

                                % Psi (basis function at quadrature points,
                                % size N X 3). Adding the scaling with
                                % g_tau and 1/(2 * volsi)
                                Psi = g_tau * (Ei * diP)' /(2 * volsi);

                                % Manually evaluating nxPsi
                                nxPsi1 = nrmx(2) * Psi(:,3) - nrmx(3) * Psi(:,2);
                                nxPsi2 = nrmx(3) * Psi(:,1) - nrmx(1) * Psi(:,3);
                                nxPsi3 = nrmx(1) * Psi(:,2) - nrmx(2) * Psi(:,1);

                              for jj = 1:Qtr % Loop over trial basis fns


                                dPsi1 = rsf_tr{jj}{1};
                                dPsi2 = rsf_tr{jj}{2};

                                dP = [dPsi1; dPsi2]; % Gradient reference element

                                dxP = DCV * dP;

                                PY1 = dxP(1) .* g_t;
                                PY2 = dxP(2) .* g_t;
                                PY3 = dxP(3) .* g_t;

                                local_matrix(ii,jj) = dot(Wh,Ker.*(nxPsi1 * PY1 + nxPsi2 * PY2 + nxPsi3 * PY3));
                              end

                          end
                          M(dofs_i(perm_i), ABC_eltj) = M(dofs_i(perm_i), ABC_eltj) + local_matrix;
                      
                      case 'nxgrad[psi]'  % Trial: nxgrad(P1), Test: nxgrad(P1)

                          for ii = 1:Qts

                            dPsi1 = rsf_ts{ii}{1};
                            dPsi2 = rsf_ts{ii}{2};

                            dP = [dPsi1; dPsi2]; % Gradient reference element

                            dxP = DCVi * dP;
                            Psix1 = dxP(1) .* g_tau;
                            Psix2 = dxP(2) .* g_tau;
                            Psix3 = dxP(3) .* g_tau;
                            
                            PX1 = nrmx(2) * Psix3 - nrmx(3) * Psix2;
                            PX2 = nrmx(3) * Psix1 - nrmx(1) * Psix3;
                            PX3 = nrmx(1) * Psix2 - nrmx(2) * Psix1;

                              for jj = 1:Qtr


                                dPsi1 = rsf_tr{jj}{1};
                                dPsi2 = rsf_tr{jj}{2};

                                dP = [dPsi1; dPsi2]; % Gradient reference element

                                dxP = DCV * dP;
                                Psiy1 = dxP(1) .* g_t;
                                Psiy2 = dxP(2) .* g_t;
                                Psiy3 = dxP(3) .* g_t;

                                PY1 = nrm(2) * Psiy3 - nrm(3) * Psiy2;
                                PY2 = nrm(3) * Psiy1 - nrm(1) * Psiy3;
                                PY3 = nrm(1) * Psiy2 - nrm(2) * Psiy1;

                                local_matrix(ii,jj) = dot(Wh,PX1 .* Ker .* PY1) ...
                                                    + dot(Wh,PX2 .* Ker .* PY2) ...
                                                    + dot(Wh,PX3 .* Ker .* PY3);
                              end

                          end
                          M(ABC_elti, ABC_eltj) = M(ABC_elti, ABC_eltj) + local_matrix;


          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      case 'n*[psi]'

                          switch ts_opr
                              
                              case '[psi]'  % Trial: ntimes(P1), Test: P0
                                  
                                  %dKer = Ker{1} * nrm(1) + Ker{2} * nrm(2)
                                  %+ Ker{3} * nrm(3); PP edit
                                  dKer = Ker(:,1) * nrm(1) + Ker(:,2) * nrm(2) + Ker(:,3) * nrm(3);
                                  for ii = 1:Qts
                                      Psix = rsf_ts{ii}(Xh) .* g_tau;

                                      for jj = 1:Qtr
                                        Psiy = rsf_tr{jj}(Yh) .* g_t;
                                        
                                        local_matrix(ii,jj) = dot(Wh,Psix .* dKer .* Psiy);
%                                         local_matrix(ii,jj) = dot(Wh,Psix .* Ker{1} .* (Psiy * nrm(1))) ...
%                                                             + dot(Wh,Psix .* Ker{2} .* (Psiy * nrm(2))) ...
%                                                             + dot(Wh,Psix .* Ker{3} .* (Psiy * nrm(3)));
                                      end

                                  end
                                  M(i, ABC_eltj) = M(i, ABC_eltj) + local_matrix;
                                  
                              case 'n*[psi]'  % Trial: ntimes(P1), Test: ntimes(P1)
                                  
                                  ndotn = nrmx(1)*nrm(1) + nrmx(2)*nrm(2) + nrmx(3)*nrm(3);
                                  for ii = 1:Qts
                                      Psix = rsf_ts{ii}(Xh) .* g_tau;

                                      for jj = 1:Qtr
                                        Psiy = rsf_tr{jj}(Yh) .* g_t;
                                        
                                        local_matrix(ii,jj) = dot(Wh, Psix .* Ker .* Psiy) * ndotn;
%                                         local_matrix(ii,jj) = dot(Wh,(Psix * nrmx(1)) .* Ker .* (Psiy * nrm(1))) ...
%                                                             + dot(Wh,(Psix * nrmx(2)) .* Ker .* (Psiy * nrm(2))) ...
%                                                             + dot(Wh,(Psix * nrmx(3)) .* Ker .* (Psiy * nrm(3)));
                                      end

                                  end
                                  M(ABC_elti, ABC_eltj) = M(ABC_elti, ABC_eltj) + local_matrix;
                                  
                                  
                                  
                          end
                          

                  
                  
                  
                  end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%%% RWG FINITE ELEMENTS BASIS FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              case 'RWG'
                  
                  
                  switch tr_opr
                      
                      case '[psi]'
                  volsi = vols(i);
                  volsj = vols(j);
                  
                  elti  = mesh.elt(i, :);
                  eltj  = mesh.elt(j, :);
                  for ii = 1:Qts
                
                
                ip1 = mod(perm_i(ii),3)+1;
                ip2 = mod(ip1,3)+1;
                
                % Flux through the face
                flux = (2*(elti(ip1) < elti(ip2))-1);
                 
                    dPsi1i = rsf_ts{ii}{1}(Xh);
                    dPsi2i = rsf_ts{ii}{2}(Xh);
                    
                    diP  = flux * [dPsi1i dPsi2i]'; % RT0 reference element
                    
                    diP = g_tau * (Ei * diP)' / (2 * volsi);
                     

                      for jj = 1:Qtr

                        jp1 = mod(perm_j(jj),3)+1;
                        jp2 = mod(jp1,3)+1;

                        % Flux through the face
                        fluxj = (2*(eltj(jp1) < eltj(jp2))-1);

                        
                        dPsi1j = rsf_tr{jj}{1}(Yh);
                        dPsi2j = rsf_tr{jj}{2}(Yh);

                        djP  = fluxj * [dPsi1j dPsi2j]'; % RT0 reference element
                        
                        djP = g_t * (E * djP)' / (2 * volsj);

                        
%                         local_matrix(ii,jj) = dot(Wh,Psi1i .* Ker .* Psi1j + Psi2i .* Ker .* Psi2j + Psi3i .* Ker .* Psi3j);
                        if size(Ker,2) == 1
                                    local_matrix(ii,jj) = dot(Wh,diP(:, 1) .* Ker .* djP(:, 1) + diP(:, 2) .* Ker .* djP(:, 2) + diP(:, 3) .* Ker .* djP(:, 3));
                        elseif size(Ker,2) ==3
                                    % Need to perform {kernel x trial}.test
                                    kerxtrial = cross(Ker,djP,2);
                                    local_matrix(ii,jj) = dot(Wh,dot(kerxtrial,diP,2));
                        end

                      end

                  end
                  % Check local2global map
                  M(dofs_i(perm_i), dofs_j(perm_j)) = M(dofs_i(perm_i), dofs_j(perm_j)) + ...
                                                      local_matrix;
                      
                      
                  end
                      
          
          
          
          end
          
          
          
%           
%           for ii = 1:Qts
%               for jj = 1:Qtr
%                   local_matrix(ii,jj) = sstri_integrate(kernel,rsf_ts{ii},rsf_tr{jj},chi_tau,chi_t,g_tau,g_t,relation);
%                   
%                   % Hard coding local to global map
%                   switch ts_typ
%                     case 'P0'
%                         II = i;
%                     case 'P1'
%                         II = ABC_elti(ii);
%                   end
%                   
%                   switch tr_typ
%                     case 'P0'
%                         JJ = j;
%                     case 'P1'
%                         JJ = ABC_eltj(jj);
%                   end
%                   
%                   M(II,JJ) = M(II,JJ) + local_matrix(ii,jj);
%                   
%               end
%           end
       
          % what's the convention here? hbasisx from which space?
          %local_integral = sstri_integrate(kernel,rsf_tr,rsf_ts,chi_tau,chi_t,g_tau,g_t,relation);
       
%        end
    end
    
    % Creating the final matrix M

end
