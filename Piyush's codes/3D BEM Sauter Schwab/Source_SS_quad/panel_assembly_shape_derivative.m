% Panel oriented assembly for double integral with a kernel in 3D

% I: index subset of elements
% J: index subset of elements
% size(I, 1) must be equal to size(J, 1).

% Pairs of elements I x J share a vertex, edge or are identical.


function M = panel_assembly_shape_derivative(mesh,kernel,trial_space,test_space, I, J,Vel,DVel)

    [X, W] = quad4D(5); Xss{1} = X{1}; Wss{1} = W{1};
    [X, W] = quad4D(5); Xss{2} = X{2}; Wss{2} = W{2};
    [X, W] = quad4D(5); Xss{3} = X{3}; Wss{3} = W{3};
    [X, W] = quad4D(5); Xss{4} = X{4}; Wss{4} = W{4};
    % Number of elements in the mesh
    N = size(mesh.elt,1);
    
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
    
    % Loop over pair of mesh elements specified in I and J
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

          
        % Evalulating the Kernel at quadrature points
        if iscell(kernel)
            Ker{1} = kernel{1}(chi_tau(Xh')',chi_t(Yh')',chi_t(Yh')'-chi_tau(Xh')');
            Ker{2} = kernel{2}(chi_tau(Xh')',chi_t(Yh')',chi_t(Yh')'-chi_tau(Xh')');
            Ker{3} = kernel{3}(chi_tau(Xh')',chi_t(Yh')',chi_t(Yh')'-chi_tau(Xh')');
        else
            Ker = kernel(chi_tau(Xh')',chi_t(Yh')',chi_t(Yh')'-chi_tau(Xh')');
        end
          

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%% Evaluating Basis Functions %%%%%%%%%%%%% 
          %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                      
                      % Hijacking the whole case for shape derivative
                      % computation (Dvel curlTg) for trial function
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
                                
                                % Components for nx grad b, grad b is
                                % constant
                                PY1 = nrm(2) * Psiy3 - nrm(3) * Psiy2;
                                PY2 = nrm(3) * Psiy1 - nrm(1) * Psiy3;
                                PY3 = nrm(1) * Psiy2 - nrm(2) * Psiy1;

                                % Evaluating DVel at the quadrature points
                                Yh_panel = chi_t(Yh'); % Size 3 X N
                                % DVi gives a matrix of size N X 3
                                % evaluating the ith row of DVel at the N
                                % quadrature points
                                DV1 = DVel{1}(Yh_panel');
                                DV2 = DVel{2}(Yh_panel');
                                DV3 = DVel{3}(Yh_panel');
                            
                                % Manually writing the matrix vector
                                % product for DVel curlTg
                                DVelcurlTg1 = sum(DV1 .* [PY1 PY2 PY3],2);
                                DVelcurlTg2 = sum(DV2 .* [PY1 PY2 PY3],2);
                                DVelcurlTg3 = sum(DV3 .* [PY1 PY2 PY3],2);

                                local_matrix(ii,jj) = dot(Wh,PX1 .* Ker .* DVelcurlTg1) ...
                                                    + dot(Wh,PX2 .* Ker .* DVelcurlTg2) ...
                                                    + dot(Wh,PX3 .* Ker .* DVelcurlTg3);
                              end

                          end
                          M(ABC_elti, ABC_eltj) = M(ABC_elti, ABC_eltj) + local_matrix;


          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      case 'n*[psi]'

                          switch ts_opr
                              
                              case '[psi]'  % Trial: ntimes(P1), Test: P0
                                  
                                  dKer = Ker{1} * nrm(1) + Ker{2} * nrm(2) + Ker{3} * nrm(3);
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
                  %disp('Case RWG');
                  
                  
                switch tr_opr
                
                case '[psi]'
                    volsi = vols(i);
                    volsj = vols(j);
                    
                    elti  = mesh.elt(i, :);
                    eltj  = mesh.elt(j, :);

                    for ii = 1:Qts % Looping over test functions

                        ip1 = mod(perm_i(ii),3)+1;
                        ip2 = mod(ip1,3)+1;
                        
                        % Flux through the face
                        flux = (2*(elti(ip1) < elti(ip2))-1);
                        
                        dPsi1i = rsf_ts{ii}{1}(Xh);
                        dPsi2i = rsf_ts{ii}{2}(Xh);
                        
                        diP  = flux * [dPsi1i dPsi2i]'; % RT0 reference element
                        
                        diP = g_tau * (Ei * diP)' / (2 * volsi);
                        
                        
                        for jj = 1:Qtr % Looping over trial functions
                        
                            jp1 = mod(perm_j(jj),3)+1;
                            jp2 = mod(jp1,3)+1;
                            
                            % Flux through the face
                            fluxj = (2*(eltj(jp1) < eltj(jp2))-1);
                            
                            
                            dPsi1j = rsf_tr{jj}{1}(Yh);
                            dPsi2j = rsf_tr{jj}{2}(Yh);
                            
                            djP  = fluxj * [dPsi1j dPsi2j]'; % RT0 reference element
                            
                            djP = g_t * (E * djP)' / (2 * volsj);
                            
                            local_matrix(ii,jj) = dot(Wh,diP(:, 1) .* Ker .* djP(:, 1) + diP(:, 2) .* Ker .* djP(:, 2) + diP(:, 3) .* Ker .* djP(:, 3));
                        
                        end
                    
                    end
                    % Check local2global map
                    M(dofs_i(perm_i), dofs_j(perm_j)) = M(dofs_i(perm_i), dofs_j(perm_j)) + ...
                                      local_matrix;

                    case 'Dvel[psi]'
                        volsi = vols(i);
                        volsj = vols(j);
                        
                        elti  = mesh.elt(i, :);
                        eltj  = mesh.elt(j, :);
    
                        for ii = 1:Qts % Looping over test functions
                            % RSF for function associated to vertex at
                            % position perm_i(ii) in elti
                            ip1 = mod(perm_i(ii),3)+1;
                            ip2 = mod(ip1,3)+1;
                            
                            % Flux through the face
                            flux = (2*(elti(ip1) < elti(ip2))-1);
                            
                            % N X 1 output vectors
                            % Xh is N X 2 input vector
                            dPsi1i = rsf_ts{ii}{1}(Xh);
                            dPsi2i = rsf_ts{ii}{2}(Xh);
                            
                            % 2 X N matrix
                            diP  = flux * [dPsi1i dPsi2i]'; % RT0 reference element
                            
                            % N X 3 matrix of the basis function evaluated
                            % at N quadrature points
                            diP = g_tau * (Ei * diP)' / (2 * volsi);
                            
                            
                            for jj = 1:Qtr % Looping over trial functions
                            
                                jp1 = mod(perm_j(jj),3)+1;
                                jp2 = mod(jp1,3)+1;
                                
                                % Flux through the face
                                fluxj = (2*(eltj(jp1) < eltj(jp2))-1);
                                
                                
                                dPsi1j = rsf_tr{jj}{1}(Yh);
                                dPsi2j = rsf_tr{jj}{2}(Yh);
                                
                                djP  = fluxj * [dPsi1j dPsi2j]'; % RT0 reference element
                                
                                % N X 3 matrix
                                djP = g_t * (E * djP)' / (2 * volsj);

                                % Evaluating DVel at the quadrature points
                                Yh_panel = chi_t(Yh'); % Size 3 X N
                                % DVi gives a matrix of size N X 3
                                % evaluating the ith row of DVel at the N
                                % quadrature points
                                DV1 = DVel{1}(Yh_panel');
                                DV2 = DVel{2}(Yh_panel');
                                DV3 = DVel{3}(Yh_panel');
                            
                                % Manually writing the matrix vector
                                % product
                                djP = [dot(DV1,djP,2) dot(DV2,djP,2) dot(DV3,djP,2)];
                                
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
          end % End switch over tr_typ
    end % End loop over panel pairs
end % Function end
