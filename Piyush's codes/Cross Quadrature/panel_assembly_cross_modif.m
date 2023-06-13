% Panel oriented assembly for double integral with a kernel in 3D

% I: index subset of elements
% J: index subset of elements
% size(I, 1) must be equal to size(J, 1).

% Pairs of elements I x J share a vertex, edge or are identical.

% trace = 0 (no trace); 1 (E x n); 2 (n x (E x n))

% mesh_test - bdry mesh
% mesh - vol mesh

function M = panel_assembly_cross_modif(mesh_test, mesh, kernel,trial_space,test_space, I, J, mult)

    if nargin < 8
       multiplier = 0;
       
    else
        multiplier = 1;
    end


    [Xss, Wss] = quad5D(5);
    
     M = zeros(test_space.ndof,trial_space.ndof);
%     M = sparse(test_space.ndof, trial_space.ndof);
%      M = spalloc(test_space.ndof, trial_space.ndof, ...
%                  floor(test_space.ndof * trial_space.ndof / 4));
   
    % Vector storing volume of the mesh elements
    vols = mesh.ndv;
    
    vols_bnd = mesh_test.ndv;
    
    normals = mesh_test.nrm;
    
    [faces, ~] = mesh.fce; fce2vtx = faces.elt;clear faces;
    [edges, ~] = mesh.edg; edg2vtx = edges.elt;clear edges;

    tr_typ = trial_space.typ;
    ts_typ = test_space.typ;

    tr_opr = trial_space.opr;
    ts_opr = test_space.opr;
    
    rsf_tr = trial_space.rsf;
    rsf_ts = test_space.rsf;
    

    [~,elt2dof_tr] = trial_space.dof;
    [~,elt2dof_ts] = test_space.dof;
    Qtr = size(elt2dof_tr,2);
    Qts = size(elt2dof_ts,2);
    
    % Convention: panel i <-> chi_tau //////<-> hbasisx
    %             panel j <-> chi_t   //////<-> hbasisy
    
    % Double loop over the mesh elements
    L = size(I, 1);
    
    
    for elt = 1:L
        
        i = I(elt);
        j = J(elt);
        

        dofs_i = elt2dof_ts(i, :);
        dofs_j = elt2dof_tr(j, :);
        measure_i = vols_bnd(i);
        normals_i = normals(i, :);
        nrmx = normals_i;
        
        
        % The Jacobian of transformation from the element in 3D to
        % reference triangle with area 1/2
 

        g_tau = 2 * measure_i;
        measure_j = vols(j);
        g_t = 6 * measure_j;

          
          % Finding the relation between panels i and j 
          [intersection, IA, IB] = intersect(mesh_test.vtx(mesh_test.elt(i, :),:),...
                                   mesh.vtx(mesh.elt(j, :),:), ...
                                   'rows', 'stable');
%           vtcs = mesh.vtx(intersection,:);
          vtcs = intersection;
          l = length(IA);
          switch l
              case 0
                  relation = "far_away";
                  % Vertices for elt i
                  Ai = mesh_test.vtx(mesh_test.elt(i,1),:)';
                  Bi = mesh_test.vtx(mesh_test.elt(i,2),:)';
                  Ci = mesh_test.vtx(mesh_test.elt(i,3),:)';
                  % Vertices for elt j
                  Aj = mesh.vtx(mesh.elt(j,1),:)';
                  Bj = mesh.vtx(mesh.elt(j,2),:)';
                  Cj = mesh.vtx(mesh.elt(j,3),:)';
                  Dj = mesh.vtx(mesh.elt(j,4),:)';
                  % Parameterizations
                  chi_tau = @(xhat) Ai + [Bi-Ai Ci-Ai]*xhat;
                  chi_t   = @(yhat) Aj + [Bj-Aj Cj-Aj Dj-Aj]*yhat;
                  ABC_elti = mesh_test.elt(i,:);
                  ABCD_eltj = mesh.elt(j,:);
                  
                  
                  perm_i = 1:3;
                  
                  perm_j = 1:4;
                         
                  miss = [1 2;1 3;1 4;2 3;2 4;3 4];
                  [~, ~, perm_ej] = intersect(sort(perm_j(miss), 2),miss, 'rows', 'stable');
                 
                  
                   
                  E = [Bj-Aj Cj-Aj Dj-Aj];
                  Etm1 = inv(E)';
                  
                  Piolaj = E / det(E);
                  
                  Ei = [Bi-Ai Ci-Ai]; EitEi = Ei' * Ei; Dxy = inv(EitEi);
                  
                  DCVx = Dxy(1, 1) * Ei(:, 1) + Dxy(1, 2) * Ei(:, 2);
                  DCVy = Dxy(2, 1) * Ei(:, 1) + Dxy(2, 2) * Ei(:, 2);
                  
                  DCVi = [DCVx DCVy];
                  
                  
                  
                  Xh = Xss{4}(:, 1:2);     
                  Yh = Xss{4}(:, 3:5);
                  Wh = Wss{4};
                 
                  
              case 1
                  relation = "common_vertex";
                  Ai = vtcs(1,:)'; Aj = Ai;
                  
                  
                  bci = mesh_test.elt(i, setdiff(1:3, IA));
                  BCi = mesh_test.vtx(bci, :);
                  Bi = BCi(1,:)';
                  Ci = BCi(2,:)';
                  
                  bcdj = mesh.elt(j, setdiff(1:4, IB));
                  BCDj = mesh.vtx(bcdj,:);
                  Bj = BCDj(1,:)';
                  Cj = BCDj(2,:)';
                  Dj = BCDj(3,:)';
                  ABC_elti = [mesh_test.elt(i, IA) bci(1) bci(2)];
                  ABCD_eltj = [mesh.elt(j, IB) bcdj(1) bcdj(2) bcdj(3)];
                  
                  
                  perm_i = [find(mesh_test.elt(i, :) == ABC_elti(1)), ...
                            find(mesh_test.elt(i, :) == ABC_elti(2)), ...
                            find(mesh_test.elt(i, :) == ABC_elti(3))];
                  
                  perm_j = [find(mesh.elt(j, :) == ABCD_eltj(1)), ...
                            find(mesh.elt(j, :) == ABCD_eltj(2)), ...
                            find(mesh.elt(j, :) == ABCD_eltj(3)), ...
                            find(mesh.elt(j, :) == ABCD_eltj(4))];
                        
                               
                  miss = [1 2;1 3;1 4;2 3;2 4;3 4];
                  [~, ~, perm_ej] = intersect(sort(perm_j(miss), 2),miss, 'rows', 'stable');
                 
                  
                  
                  % 0,0 for common point
                  chi_tau = @(xhat) Ai + [Bi-Ai Ci-Ai]*xhat;
                  chi_t = @(yhat) Aj + [Bj-Aj Cj-Aj Dj-Aj]*yhat;
                  
                   
                  E = [Bj-Aj Cj-Aj Dj-Aj];
                  Etm1 = inv(E)';
                  
                  Piolaj = E / det(E);
                  
                  Ei = [Bi-Ai Ci-Ai]; EitEi = Ei' * Ei; Dxy = inv(EitEi);
                  
                  DCVx = Dxy(1, 1) * Ei(:, 1) + Dxy(1, 2) * Ei(:, 2);
                  DCVy = Dxy(2, 1) * Ei(:, 1) + Dxy(2, 2) * Ei(:, 2);
                  
                  DCVi = [DCVx DCVy];
                  
                 
                  
                  Xh = Xss{3}(:, 1:2);
                  Yh = Xss{3}(:, 3:5);
                  
%                   
%                   Xh(:,1) = Xh(:,1) - Xh(:,2);
%                   Yh(:,1) = Yh(:,1) - Yh(:,2) - Yh(:,3);
                  
                  Wh = Wss{3};
                  
              case 2
                  relation = "common_edge";
                  Ai = vtcs(1,:)'; Aj = Ai;
                  Bi = vtcs(2,:)'; Bj = Bi;
                  
                  
                  
                  ci = mesh_test.elt(i, setdiff(1:3, IA));
                  Ci = mesh_test.vtx(ci, :)';
                  
                  
                  cdj = mesh.elt(j, setdiff(1:4, IB));
                  
                  Cj = mesh.vtx(cdj(1), :)';
                  Dj = mesh.vtx(cdj(2), :)';
                  ABC_elti = [mesh_test.elt(i, IA) ci];
                  ABCD_eltj = [mesh.elt(j, IB) cdj(1) cdj(2)];
                  
                  
                  perm_i = [find(mesh_test.elt(i, :) == ABC_elti(1)), ...
                            find(mesh_test.elt(i, :) == ABC_elti(2)), ...
                            find(mesh_test.elt(i, :) == ABC_elti(3))];
                  
                  perm_j = [find(mesh.elt(j, :) == ABCD_eltj(1)), ...
                            find(mesh.elt(j, :) == ABCD_eltj(2)), ...
                            find(mesh.elt(j, :) == ABCD_eltj(3)), ...
                            find(mesh.elt(j, :) == ABCD_eltj(4))];
                        
                               
                  miss = [1 2;1 3;1 4;2 3;2 4;3 4];
                  [~, ~, perm_ej] = intersect(sort(perm_j(miss), 2),miss, 'rows', 'stable');
                 
                  
                  
                  % 0,0 and 1,0 for common points
                  chi_tau = @(xhat) Ai + [Bi-Ai Ci-Ai]*xhat;
                  chi_t   = @(yhat) Aj + [Bj-Aj Cj-Aj Dj-Aj]*yhat;
                  
                   
                  E = [Bj-Aj Cj-Aj Dj-Aj];
                  Etm1 = inv(E)';
                  
                  Piolaj = E / det(E);
                  
                  Ei = [Bi-Ai Ci-Ai]; EitEi = Ei' * Ei; Dxy = inv(EitEi);
                  
                  DCVx = Dxy(1, 1) * Ei(:, 1) + Dxy(1, 2) * Ei(:, 2);
                  DCVy = Dxy(2, 1) * Ei(:, 1) + Dxy(2, 2) * Ei(:, 2);
                  
                  DCVi = [DCVx DCVy];
                  
                 
                  
                  
                  Xh = Xss{2}(:, 1:2);              
                  Yh = Xss{2}(:, 3:5);
%                   
%                   Xh(:,1) = Xh(:,1) - Xh(:,2);
%                   Yh(:,1) = Yh(:,1) - Yh(:,2) - Yh(:,3);
                  Wh = Wss{2};
                  
              case 3
                  relation = "common_face";
                  Aj = vtcs(1,:)'; Ai = Aj;
                  Bj = vtcs(2,:)'; Bi = Bj;
                  Cj = vtcs(3,:)'; Ci = Cj;
                  
                  dj = mesh.elt(j, setdiff(1:4, IB));
                  Dj = mesh.vtx(dj, :)';
                  chi_tau = @(xhat) Ai + [Bi-Ai Ci-Ai]*xhat;
                  chi_t = @(yhat) Aj + [Bj-Aj Cj-Aj Dj-Aj]*yhat;
                  ABC_elti = [mesh_test.elt(i, IA)];
                  
                  ABCD_eltj = [mesh.elt(j, IB) dj];
                  
                  
                  perm_i = [find(mesh_test.elt(i, :) == ABC_elti(1)), ...
                            find(mesh_test.elt(i, :) == ABC_elti(2)), ...
                            find(mesh_test.elt(i, :) == ABC_elti(3))];
                  
                  perm_j = [find(mesh.elt(j, :) == ABCD_eltj(1)), ...
                            find(mesh.elt(j, :) == ABCD_eltj(2)), ...
                            find(mesh.elt(j, :) == ABCD_eltj(3)), ...
                            find(mesh.elt(j, :) == ABCD_eltj(4))];
                        
                         
                  miss = [1 2;1 3;1 4;2 3;2 4;3 4];
                  [~, ~, perm_ej] = intersect(sort(perm_j(miss), 2),miss, 'rows', 'stable');
                 
                  
                  
                   
                  E = [Bj-Aj Cj-Aj Dj-Aj];
                  Etm1 = inv(E)';
                  
                  Piolaj = E / det(E);
                 
                  Ei = [Bi-Ai Ci-Ai]; EitEi = Ei' * Ei; Dxy = inv(EitEi);
                  
                  DCVx = Dxy(1, 1) * Ei(:, 1) + Dxy(1, 2) * Ei(:, 2);
                  DCVy = Dxy(2, 1) * Ei(:, 1) + Dxy(2, 2) * Ei(:, 2);
                  
                  DCVi = [DCVx DCVy];
                  
                  
                  
                  Xh = Xss{1}(:, 1:2);     
                  Yh = Xss{1}(:, 3:5);
                  
%                   
%                   Xh(:,1) = Xh(:,1) - Xh(:,2);
%                   Yh(:,1) = Yh(:,1) - Yh(:,2) - Yh(:,3);
                  Wh = Wss{1};
              
          end
          
          
          local_matrix = zeros(Qts,Qtr);
          
%           g_tau = sqrt(det(Ei' * Ei));
          
          

          if iscell(kernel)
            Ker{1} = kernel{1}(chi_tau(Xh')',chi_t(Yh')',chi_t(Yh')'-chi_tau(Xh')');
            Ker{2} = kernel{2}(chi_tau(Xh')',chi_t(Yh')',chi_t(Yh')'-chi_tau(Xh')');
            Ker{3} = kernel{3}(chi_tau(Xh')',chi_t(Yh')',chi_t(Yh')'-chi_tau(Xh')');
          
          else
              
            Ker = kernel(chi_tau(Xh')',chi_t(Yh')',chi_t(Yh')'-chi_tau(Xh')');
          end       
          
          
          if multiplier
             
              if iscell(mult)
                 mult1 = mult{1}(chi_t(Yh')');
                 mult2 = mult{2}(chi_t(Yh')');
                 mult3 = mult{3}(chi_t(Yh')');
                  
              else
                 Mult = mult(chi_t(Yh')');
              end
              
          end
          
          
          %---------------------------------------------------------------%
          %---------------------------------------------------------------%
          switch tr_typ
              case 'P0'
                  switch ts_typ
                      case 'P0'
                          assert(strcmp(ts_opr,'n*[psi]')); % Vectorial kernel is hard coded in this case
                          Kerdotn = sum(normals_i .* Ker,2);
                          for ii = 1:Qts
                              Psix = rsf_ts{ii}(Xh) .* g_tau;

                              for jj = 1:Qtr
                                Psiy = rsf_tr{jj}(Yh) .* g_t;
                                
                                local_matrix(ii,jj) = dot(Wh, Psix .* Kerdotn .* Psiy);
                              end

                          end
                          M(ABC_elti, ABCD_eltj) = M(ABC_elti, ABCD_eltj) + local_matrix;
                          % Alternative
                          %M(i, j) = M(i, j) + local_matrix;

                      case 'P1'
                          assert(strcmp(ts_opr,'[psi]')); 
                          for ii = 1:Qts
                              Psix = rsf_ts{ii}(Xh) .* g_tau;

                              for jj = 1:Qtr
                                Psiy = rsf_tr{jj}(Yh) .* g_t;
                                
                                local_matrix(ii,jj) = dot(Wh, Psix .* Kerdotn .* Psiy);
                              end

                          end
                          M(ABC_elti, ABCD_eltj) = M(ABC_elti, ABCD_eltj) + local_matrix;
                          % Alternative
                          M(ABC_elti, j) = M(ABC_elti, j) + local_matrix;
                  end

              case 'P1'
              
          %---------------------------------------------------------------%
          %---------------------------------------------------------------%
          switch tr_opr
              case '[psi]'
                  
                  switch ts_opr
                      case '[psi]'
                          
          %---------------------------------------------------------------%
                          for ii = 1:Qts
                              Psix = rsf_ts{ii}(Xh) .* g_tau;

                              for jj = 1:Qtr
                                Psiy = rsf_tr{jj}(Yh) .* g_t;
                                local_matrix(ii,jj) = dot(Wh,Psix .* Ker .* Psiy);
                              end

                          end
                          M(i, ABCD_eltj) = M(i, ABCD_eltj) + local_matrix;

          %---------------------------------------------------------------%
                      case 'n*[psi]'   

                          dKer = normals_i(1) * Ker{1} + normals_i(2) * Ker{2} + normals_i(3) * Ker{3}; 
                          for ii = 1:Qts
                              Psix = rsf_ts{ii}(Xh) .* g_tau;

                              for jj = 1:Qtr
                                Psiy = rsf_tr{jj}(Yh) .* g_t;
                                
                                local_matrix(ii,jj) = dot(Wh, Psix .* dKer .* Psiy);
%                                 local_matrix(ii,jj) = dot(Wh,normals_i(1) .* Psix .* Ker{1} .* Psiy) + ...
%                                                       dot(Wh,normals_i(2) .* Psix .* Ker{2} .* Psiy) + ...
%                                                       dot(Wh,normals_i(3) .* Psix .* Ker{3} .* Psiy);
                              end

                          end
                          % Hard-coded for P1 volume, P0 boundary.
                          M(ABC_elti, ABCD_eltj) = M(ABC_elti, ABCD_eltj) + local_matrix;
                          
          %---------------------------------------------------------------%
                      case 'nxgrad[psi]'   
%                           
%                           nxKer{1} = normals_i(2) * Ker{3} - normals_i(3) * Ker{2};
%                           nxKer{2} = normals_i(3) * Ker{1} - normals_i(1) * Ker{3};
%                           nxKer{3} = normals_i(1) * Ker{2} - normals_i(2) * Ker{1};

                          for ii = 1:Qts
                            dPsi1 = rsf_ts{ii}{1} * g_tau;
                            dPsi2 = rsf_ts{ii}{2} * g_tau;


                            dP = [dPsi1; dPsi2]; % Gradient reference element

                            dxP = DCVi * dP;

                            Psix1 = dxP(1);
                            Psix2 = dxP(2);
                            Psix3 = dxP(3);
                            
                             
                            PX1 = normals_i(2) * Psix3 - normals_i(3) * Psix2;
                            PX2 = normals_i(3) * Psix1 - normals_i(1) * Psix3;
                            PX3 = normals_i(1) * Psix2 - normals_i(2) * Psix1;


                              for jj = 1:Qtr
                                Psiy = rsf_tr{jj}(Yh) .* g_t;
%                                 local_matrix(ii,jj) = dot(Wh,PX1 .* nxKer{1} .* Psiy) + ...
%                                                       dot(Wh,PX2 .* nxKer{2} .* Psiy) + ...
%                                                       dot(Wh,PX3 .* nxKer{3} .* Psiy);

                                local_matrix(ii,jj) = dot(Wh,(PX1 .* Ker{1} + PX2 .* Ker{2} + PX3 .* Ker{3}) .* Psiy);
                              end

                          end
                          
                          M(ABC_elti, ABCD_eltj) = M(ABC_elti, ABCD_eltj) + local_matrix;

                  end
                  
          %---------------------------------------------------------------%
          %---------------------------------------------------------------%
              case 'grad[psi]'
                  
                  
                  for ii = 1:Qts
                      Psix = rsf_ts{ii}(Xh) .* g_tau;

                      for jj = 1:Qtr
                        dPsi1 = rsf_tr{jj}{1}(Yh);
                        dPsi2 = rsf_tr{jj}{2}(Yh);
                        dPsi3 = rsf_tr{jj}{3}(Yh);

                        dP = [dPsi1; dPsi2; dPsi3]; % Gradient reference element

                        dxP = Etm1 * dP;
                        Psiy1 = dxP(1) .* g_t;
                        Psiy2 = dxP(2) .* g_t;
                        Psiy3 = dxP(3) .* g_t;
                        
                        
                        local_matrix(ii,jj) = dot(Wh,Psix .* Ker{1} .* Psiy1) ...
                                            + dot(Wh,Psix .* Ker{2} .* Psiy2) ...
                                            + dot(Wh,Psix .* Ker{3} .* Psiy3);
                      end

                  end
                  M(ABCD_elti, ABCD_eltj) = M(ABCD_elti, ABCD_eltj) + ...
                                            local_matrix;

          end
          
              case 'RWG'
                  switch tr_opr
                      case '[psi]'
                          for ii = 1:Qts
                ip1 = mod(ii,4)+1;
                ip2 = mod(ip1,4)+1;
                ip3 = mod(ip2,4)+1;
                
                % Flux through the face
                flux = (2*(mesh.elt(i,perm_i(ip1)) < mesh.elt(i,perm_i(ip2)))-1) ...
                    .*(2*(mesh.elt(i,perm_i(ip2)) < mesh.elt(i,perm_i(ip3)))-1) ...
                    .*(2*(mesh.elt(i,perm_i(ip1)) < mesh.elt(i,perm_i(ip3)))-1) ...
                    ./0.5;
                
%                 flux_old = (2*(mesh.elt(i,ip1) < mesh.elt(i,ip2))-1) ...
%                     .*(2*(mesh.elt(i,ip2) < mesh.elt(i,ip3))-1) ...
%                     .*(2*(mesh.elt(i,ip1) < mesh.elt(i,ip3))-1) ...
%                     ./0.5;
                 
                
                if mod(ii,2)==0
                    flux = - flux;
                end
                
                      
                    dPsi1i = rsf_tr{ii}{1}(Yh);
                    dPsi2i = rsf_tr{ii}{2}(Yh);
                    dPsi3i = rsf_tr{ii}{3}(Yh);
                    
                    diP  = flux * [dPsi1i dPsi2i dPsi3i]'; % RT0 reference element
                    diP = (Piolai * diP)';
                    
                    Psi1i = diP(:, 1) .* g_tau;
                    Psi2i = diP(:, 2) .* g_tau;
                    Psi3i = diP(:, 3) .* g_tau;

                      for jj = 1:Qtr
                          
                        jp1 = mod(jj,4)+1;
                        jp2 = mod(jp1,4)+1;
                        jp3 = mod(jp2,4)+1;

                        % Flux through the face
                        fluxj = (2*(mesh.elt(j,perm_j(jp1)) < mesh.elt(j,perm_j(jp2)))-1) ...
                            .*(2*(mesh.elt(j,perm_j(jp2)) < mesh.elt(j,perm_j(jp3)))-1) ...
                            .*(2*(mesh.elt(j,perm_j(jp1)) < mesh.elt(j,perm_j(jp3)))-1) ...
                            ./0.5;
%               

                        if mod(jj,2)==0
                            fluxj = - fluxj;
                        end
                        dPsi1j = rsf_tr{jj}{1}(Yh);
                        dPsi2j = rsf_tr{jj}{2}(Yh);
                        dPsi3j = rsf_tr{jj}{3}(Yh);

                        djP  = fluxj * [dPsi1j dPsi2j dPsi3j]'; % RT0 reference element
                        djP = (Piolaj * djP)';

                        Psi1j = djP(:,1) .* g_t;
                        Psi2j = djP(:,2) .* g_t;
                        Psi3j = djP(:,3) .* g_t;
                        
                        
                        local_matrix(ii,jj) = dot(Wh,Psi1i .* Ker .* Psi1j) ...
                                            + dot(Wh,Psi2i .* Ker .* Psi2j) ...
                                            + dot(Wh,Psi3i .* Ker .* Psi3j);
                      end

                  end
                  
                  M(dofs_i(perm_i), dofs_j(perm_j)) = M(dofs_i(perm_i), dofs_j(perm_j)) + ...
                                                      local_matrix;

                      case 'div[psi]'

                          for ii = 1:Qts
                              Psix = rsf_ts{ii}(Xh) .* g_tau;

                              for jj = 1:Qtr
                                  jp1 = mod(jj,4)+1;
                                  jp2 = mod(jp1,4)+1;
                                  jp3 = mod(jp2,4)+1;

                                  % Flux through the face
                                  fluxj = (2*(mesh.elt(j,perm_j(jp1)) < mesh.elt(j,perm_j(jp2)))-1) ...
                                        .*(2*(mesh.elt(j,perm_j(jp2)) < mesh.elt(j,perm_j(jp3)))-1) ...
                                        .*(2*(mesh.elt(j,perm_j(jp1)) < mesh.elt(j,perm_j(jp3)))-1) ...
                                        ./0.5;
%               

                                    if mod(jj,2)==0
                                        fluxj = - fluxj;
                                    end
%                                     dPsi1j = rsf_tr{jj}{1}(Yh);
%                                     dPsi2j = rsf_tr{jj}{2}(Yh);
%                                     dPsi3j = rsf_tr{jj}{3}(Yh);
%             
%                                     djP  = fluxj * [dPsi1j dPsi2j dPsi3j]'; % RT0 reference element

                                    div_djP0 = fluxj * 3; % div RT0 reference element

                                    %djP = (Piolaj * djP)';

                                    % Transforming to the non reference
                                    % system
                                    div_djP = div_djP0/det(E);
            
                                    div_djP = div_djP .* g_t;

                                    if size(Ker,2) == 3 % Vectorial kernel, test space ntimes(P1) used
                                        Kerdotn = normals_i(1) * Ker(:,1) + normals_i(2) * Ker(:,2) + normals_i(3) * Ker(:,3); 
                                        local_matrix(ii,jj) = dot(Wh,div_djP * Kerdotn .*Psix);
                                    else % test space P0 assumed
                                        local_matrix(ii,jj) = dot(Wh,div_djP * Ker .*Psix);
                                    end
                              end
                          end
                          if size(Ker,2) == 3 % Vectorial kernel, test space ntimes(P1) used
                                M(ABC_elti,dofs_j(perm_j)) = M(ABC_elti,dofs_j(perm_j)) + local_matrix;
                          else % test space P0 assumed
                                M(i,dofs_j(perm_j)) = M(i,dofs_j(perm_j)) + local_matrix;
                          end

                  end
                  
                  
          
                                                  
                                                  
              case 'NED'
                  
                  switch ts_opr
                      case '[psi]'
                  
                  switch tr_opr
                      case '[psi]'
                        
                  volsi = vols_bnd(i);
                  
                  elti  = mesh_test.elt(i, :);
                  eltj  = mesh.elt(j, :);
                  
                  ej = ABCD_eltj(miss);
                  for ii = 1:Qts
                
                
                ip1 = mod(perm_i(ii),3)+1;
                ip2 = mod(ip1,3)+1;
                
                % Flux through the face
                flux = (2*(elti(ip1) < elti(ip2))-1);
                 
                    dPsi1i = rsf_ts{ii}{1}(Xh) * g_tau;
                    dPsi2i = rsf_ts{ii}{2}(Xh) * g_tau;
                    
                    diP  = flux * [dPsi1i dPsi2i]'; % RT0 reference element
                    
                    diP = (Ei * diP)' / (2 * volsi);
                     

                      for jj = 1:Qtr

                      
                          sign_j = 2 * (ej(jj, 1) < ej(jj, 2)) -1;
                          dPsi1 = rsf_tr{jj}{1}(Yh) * g_t;
                          dPsi2 = rsf_tr{jj}{2}(Yh) * g_t;
                          dPsi3 = rsf_tr{jj}{3}(Yh) * g_t;
                           
                          dP = [dPsi1 dPsi2 dPsi3]'; 
                           
                          dxP = sign_j * (Etm1 * dP)';
                          Psiy1 = dxP(:, 1);
                          Psiy2 = dxP(:, 2);
                          Psiy3 = dxP(:, 3);

                        if multiplier && ~iscell(kernel) && ~iscell(mult)   % (N(pE), beta)
                        
                            local_matrix(ii,jj) = dot(Wh,Mult .* Ker .* (diP(:, 1) .* Psiy1  + ...
                                                                      diP(:, 2) .* Psiy2  + ...
                                                                      diP(:, 3) .* Psiy3));
                            
                        elseif multiplier && ~iscell(kernel) && iscell(mult) % (N(nabla p x E), beta)
                            
                            local_matrix(ii,jj) = dot(Wh,Ker .* (diP(:, 1) .* (mult2 .* Psiy3 - mult3 .* Psiy2) + ...
                                                          diP(:, 2) .* (mult3 .* Psiy1 - mult1 .* Psiy3) + ...
                                                          diP(:, 3) .* (mult1 .* Psiy2 - mult2 .* Psiy1) ));

                                                      
                        elseif multiplier && iscell(kernel) && iscell(mult) % (grad N( grad p dot E), beta)
                            
                            local_matrix(ii,jj) = dot(Wh,(diP(:, 1) .* Ker{1}  + ...
                                                          diP(:, 2) .* Ker{2}  + ...
                                                          diP(:, 3) .* Ker{3} ).* (mult1 .* Psiy1 + mult2 .* Psiy2 + mult3 .* Psiy3) );

                        else
                            
                            local_matrix(ii,jj) = dot(Wh,Ker .* (diP(:, 1) .* Psiy1 + diP(:, 2) .* Psiy2 + diP(:, 3) .* Psiy3));
                        
                        end
                      end

                  end
                  
                  M(dofs_i(perm_i), dofs_j(perm_ej)) = full(M(dofs_i(perm_i), dofs_j(perm_ej))) + ...
                                                      local_matrix;
                                   
                 
                      
                      case 'curl[psi]' % Trial : curl(NED), test - RWG
                        assert(strcmp(ts_typ,'RWG'));
                  volsi = vols_bnd(i);
                  
                  elti  = mesh_test.elt(i, :);
                  eltj  = mesh.elt(j, :);
                  
                  ej = ABCD_eltj(miss);
                  for ii = 1:Qts
                
                
                ip1 = mod(perm_i(ii),3)+1;
                ip2 = mod(ip1,3)+1;
                
                % Flux through the face
                flux = (2*(elti(ip1) < elti(ip2))-1);
                 
                    dPsi1i = rsf_ts{ii}{1}(Xh) * g_tau;
                    dPsi2i = rsf_ts{ii}{2}(Xh) * g_tau;
                    
                    diP  = flux * [dPsi1i dPsi2i]'; % RT0 reference element
                    
                    diP = (Ei * diP)' / (2 * volsi);
                     

                      for jj = 1:Qtr

                      
                          sign_j = 2 * (ej(jj, 1) < ej(jj, 2)) -1;
                          dPsi1 = rsf_tr{jj}{1};
                          dPsi2 = rsf_tr{jj}{2};
                          dPsi3 = rsf_tr{jj}{3};
                           
                          dP = [dPsi1 dPsi2 dPsi3]'; 
                           
                          dxP = 2 * sign_j * (E * dP)' / (det(E));
                          Psiy1 = dxP(:, 1);
                          Psiy2 = dxP(:, 2);
                          Psiy3 = dxP(:, 3);

                          if size(Ker,2)==3 % (gradG X NED).RWG
                            % computing Ker X Psiy
                            KerxPsiy_1 = Ker(:,2) * Psiy3 - Ker(:,3) * Psiy2;
                            KerxPsiy_2 = Ker(:,3) * Psiy1 - Ker(:,1) * Psiy3;
                            KerxPsiy_3 = Ker(:,1) * Psiy2 - Ker(:,2) * Psiy1;
                            KerxPsi = [KerxPsiy_1 KerxPsiy_2 KerxPsiy_3];
                            local_matrix(ii,jj) = dot(Wh,dot(KerxPsi,diP,2),1);

                          else
                              local_matrix(ii,jj) = dot(Wh,diP(:, 1) .* Ker .* Psiy1 + diP(:, 2) .* Ker .* Psiy2 + diP(:, 3) .* Ker .* Psiy3);
                          end

%                         if multiplier && ~iscell(kernel) && ~iscell(mult)   % (N(p curlE), beta)
%                         
%                             local_matrix(ii,jj) = dot(Wh,Mult .* (diP(:, 1) .* Ker .* Psiy1  + ...
%                                                                       diP(:, 2) .* Ker .* Psiy2  + ...
%                                                                       diP(:, 3) .* Ker .* Psiy3));
%                             
%                         elseif multiplier && ~iscell(kernel) && iscell(mult) % (N(nabla p x curlE), beta)
%                             
%                             local_matrix(ii,jj) = dot(Wh,(diP(:, 1) .* Ker .* (mult2 .* Psiy3 - mult3 .* Psiy2) + ...
%                                                           diP(:, 2) .* Ker .* (mult3 .* Psiy1 - mult1 .* Psiy3) + ...
%                                                           diP(:, 3) .* Ker .* (mult1 .* Psiy2 - mult2 .* Psiy1) ));
% 
%                                                       
%                         elseif multiplier && iscell(kernel) && iscell(mult) % (grad N( grad p dot curlE), beta)
%                             
%                             local_matrix(ii,jj) = dot(Wh,(diP(:, 1) .* Ker{1}  + ...
%                                                           diP(:, 2) .* Ker{2}  + ...
%                                                           diP(:, 3) .* Ker{3}) .* (mult1 .* Psiy1 + mult2 .* Psiy2 + mult3 .* Psiy3) );
% 
%                         else
%                             
%                             local_matrix(ii,jj) = dot(Wh,diP(:, 1) .* Ker .* Psiy1 + diP(:, 2) .* Ker .* Psiy2 + diP(:, 3) .* Ker .* Psiy3);
%                         
%                         end
                      end

                  end
                  
                  M(dofs_i(perm_i), dofs_j(perm_ej)) = M(dofs_i(perm_i), dofs_j(perm_ej)) + ...
                                                      local_matrix;
                      
                  end
                  
                      case 'nx[psi]'
                     
                  switch tr_opr
                      case '[psi]'
                        
                  volsi = vols_bnd(i);
                  
                  elti  = mesh_test.elt(i, :);
                  eltj  = mesh.elt(j, :);
                  
                  ej = ABCD_eltj(miss);
                  for ii = 1:Qts
                
                
                ip1 = mod(perm_i(ii),3)+1;
                ip2 = mod(ip1,3)+1;
                
                % Flux through the face
                flux = (2*(elti(ip1) < elti(ip2))-1);
                 
                    dPsi1i = rsf_ts{ii}{1}(Xh) * g_tau;
                    dPsi2i = rsf_ts{ii}{2}(Xh) * g_tau;
                    
                    diP  = flux * [dPsi1i dPsi2i]'; % RT0 reference element
                    
                    diP = (Ei * diP)' / (2 * volsi);
                    

                    PX1 = nrmx(2) * diP(:, 3) - nrmx(3) * diP(:, 2);
                    PX2 = nrmx(3) * diP(:, 1) - nrmx(1) * diP(:, 3);
                    PX3 = nrmx(1) * diP(:, 2) - nrmx(2) * diP(:, 1);
                     

                      for jj = 1:Qtr

                      
                          sign_j = 2 * (ej(jj, 1) < ej(jj, 2)) -1;
                          dPsi1 = rsf_tr{jj}{1}(Yh) * g_t;
                          dPsi2 = rsf_tr{jj}{2}(Yh) * g_t;
                          dPsi3 = rsf_tr{jj}{3}(Yh) * g_t;
                           
                          dP = [dPsi1 dPsi2 dPsi3]'; 
                           
                          dxP = sign_j * (Etm1 * dP)';
                          Psiy1 = dxP(:, 1);
                          Psiy2 = dxP(:, 2);
                          Psiy3 = dxP(:, 3);

                        if multiplier && ~iscell(kernel) && ~iscell(mult)   % (N(pE), beta)
                        
                            local_matrix(ii,jj) = dot(Wh,Mult .* (PX1 .* Ker .* Psiy1  + ...
                                                                      PX2 .* Ker .* Psiy2  + ...
                                                                      PX3 .* Ker .* Psiy3));
                            
                        elseif multiplier && ~iscell(kernel) && iscell(mult) % (N(nabla p x E), beta)
                            
                            local_matrix(ii,jj) = dot(Wh,(PX1 .* Ker .* (mult2 .* Psiy3 - mult3 .* Psiy2) + ...
                                                          PX2 .* Ker .* (mult3 .* Psiy1 - mult1 .* Psiy3) + ...
                                                          PX3 .* Ker .* (mult1 .* Psiy2 - mult2 .* Psiy1) ));

                                                      
                        elseif multiplier && iscell(kernel) && iscell(mult) % (grad N( grad p dot E), beta)
                            
                            local_matrix(ii,jj) = dot(Wh,(PX1 .* Ker{1} + ...
                                                          PX2 .* Ker{2} + ...
                                                          PX3 .* Ker{3}) .* (mult1 .* Psiy1 + mult2 .* Psiy2 + mult3 .* Psiy3) );

                        else
                            
                            local_matrix(ii,jj) = dot(Wh,PX1 .* Ker .* Psiy1 + PX2 .* Ker .* Psiy2 + PX3 .* Ker .* Psiy3);
                        
                        end
                      end

                  end
                  
                  M(dofs_i(perm_i), dofs_j(perm_ej)) = M(dofs_i(perm_i), dofs_j(perm_ej)) + ...
                                                      local_matrix;
                                   
                 
                      
                      case 'curl[psi]'
                        
                  volsi = vols_bnd(i);
                  
                  elti  = mesh_test.elt(i, :);
                  eltj  = mesh.elt(j, :);
                  
                  ej = ABCD_eltj(miss);
                  for ii = 1:Qts
                
                
                ip1 = mod(perm_i(ii),3)+1;
                ip2 = mod(ip1,3)+1;
                
                % Flux through the face
                flux = (2*(elti(ip1) < elti(ip2))-1);
                 
                    dPsi1i = rsf_ts{ii}{1}(Xh) * g_tau;
                    dPsi2i = rsf_ts{ii}{2}(Xh) * g_tau;
                    
                    diP  = flux * [dPsi1i dPsi2i]'; % RT0 reference element
                    
                    diP = (Ei * diP)' / (2 * volsi);
                    
                    
                    PX1 = nrmx(2) * diP(:, 3) - nrmx(3) * diP(:, 2);
                    PX2 = nrmx(3) * diP(:, 1) - nrmx(1) * diP(:, 3);
                    PX3 = nrmx(1) * diP(:, 2) - nrmx(2) * diP(:, 1);
                     

                      for jj = 1:Qtr

                      
                          sign_j = 2 * (ej(jj, 1) < ej(jj, 2)) -1;
                          dPsi1 = rsf_tr{jj}{1};
                          dPsi2 = rsf_tr{jj}{2};
                          dPsi3 = rsf_tr{jj}{3};
                           
                          dP = [dPsi1 dPsi2 dPsi3]'; 
                           
                          dxP = 2 * sign_j * (E * dP)' / (det(E));
                          Psiy1 = dxP(:, 1);
                          Psiy2 = dxP(:, 2);
                          Psiy3 = dxP(:, 3);


                        if multiplier && ~iscell(kernel) && ~iscell(mult)   % (N(p curlE), beta)
                        
                            local_matrix(ii,jj) = dot(Wh,Mult .* (PX1 .* Ker .* Psiy1  + ...
                                                                      PX2 .* Ker .* Psiy2  + ...
                                                                      PX3 .* Ker .* Psiy3));
                            
                        elseif multiplier && ~iscell(kernel) && iscell(mult) % (N(nabla p x curlE), beta)
                            
                            local_matrix(ii,jj) = dot(Wh,(PX1 .* Ker .* (mult2 .* Psiy3 - mult3 .* Psiy2) + ...
                                                          PX2 .* Ker .* (mult3 .* Psiy1 - mult1 .* Psiy3) + ...
                                                          PX3 .* Ker .* (mult1 .* Psiy2 - mult2 .* Psiy1) ));

                                                      
                        elseif multiplier && iscell(kernel) && iscell(mult) % (grad N( grad p dot curlE), beta)
                            
                            local_matrix(ii,jj) = dot(Wh,(PX1 .* Ker{1} + ...
                                                          PX2 .* Ker{2} + ...
                                                          PX3 .* Ker{3}) .* (mult1 .* Psiy1 + mult2 .* Psiy2 + mult3 .* Psiy3) );

                        else
                            
                            local_matrix(ii,jj) = dot(Wh,PX1 .* Ker .* Psiy1 + PX2 .* Ker .* Psiy2 + PX3 .* Ker .* Psiy3);
                        
                        end
                      end

                  end
                  
                  M(dofs_i(perm_i), dofs_j(perm_ej)) = M(dofs_i(perm_i), dofs_j(perm_ej)) + ...
                                                      local_matrix;
                      
                  end
                  
                  
                          
                          
                  end
                          
                          
                      case 'div[psi]'
                  
                  
                  
          end
          
          
    end
    

end
