% Panel oriented assembly for double integral with a kernel in 3D

function M = panel_oriented_assembly(mesh,kernel,trial_space,test_space)
    % Number of elements in the mesh
    N = mesh.elt;
    
    M = zeros(test_space.ndof,trial_space.ndof);
    
    % Convention: panel i <-> chi_tau //////<-> hbasisx
    %             panel j <-> chi_t   //////<-> hbasisy
    
    % Double loop over the mesh elements
    for i = 1:N
        area_i = mesh.ndv(i);
        g_tau = @(x) area_i;
        
       for j = 1:N
          area_j = mesh.ndv(j);
          g_t = @(y) area_j;
          
          % Finding the relation between panels i and j 
          intersection = intersect(mesh.elt(i,:),mesh.elt(j,:));
          vtcs = mesh.vtx(intersection,:);
          l = length(intersection);
          switch l
              case 0
                  relation = "far";
                  Ai = mesh.elt(i,1);
                  Bi = mesh.elt(i,2);
                  Ci = mesh.elt(i,3);
                  Aj = mesh.elt(j,1);
                  Bj = mesh.elt(j,2);
                  Cj = mesh.elt(j,3);
                  chi_tau = @(xhat) Ai + [Bi-Ai Ci-Bi]*xhat;
                  chi_t = @(yhat) Aj + [Bj-Aj Cj-Bj]*yhat;
                  ABC_elti = mesh.elt(i,:);
                  ABC_eltj = mesh.elt(j,:);
                  
              case 1
                  relation = "common_vertex";
                  Ai = vtcs(1); Aj = Ai;
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
                  chi_tau = @(xhat) Ai + [Bi-Ai Ci-Bi]*xhat;
                  chi_t = @(yhat) Aj + [Bj-Aj Cj-Bj]*yhat;
                  
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
                  chi_tau = @(xhat) Ai + [Bi-Ai Ci-Bi]*xhat;
                  chi_t = @(yhat) Aj + [Bj-Aj Cj-Bj]*yhat;
                  
              case 3
                  relation = "identical";
                  Aj = vtcs(1,:)'; Ai = Aj;
                  Bj = vtcs(2,:)'; Bi = Bj;
                  Cj = vtcs(3,:)'; Ci = Cj;
                  chi_tau = @(xhat) Ai + [Bi-Ai Ci-Bi]*xhat;
                  chi_t = @(yhat) Aj + [Bj-Aj Cj-Bj]*yhat;
                  ABC_elti = intersection;
                  ABC_eltj = intersection;
          end
          
          [~,elt2dof_tr] = trial_space.dof;
          [~,elt2dof_ts] = test_space.dof;
          Qtr = size(elt2dof_tr,2);
          Qts = size(elt2dof_ts,2);
          rsf_tr = trial_space.rsf;
          rsf_ts = test_space.rsf;
          tr_typ = trial_space.typ;
          ts_typ = test_space.typ;
          
          local_matrix = zeros(Qtr,Qts);
          
          for ii = 1:Qts
              for jj = 1:Qtr
                  local_matrix(ii,jj) = sstri_integrate(kernel,rsf_ts{ii},rsf_tr{jj},chi_tau,chi_t,g_tau,g_t,relation);
                  
                  % Hard coding local to global map
                  switch ts_typ
                    case 'P0'
                        II = i;
                    case 'P1'
                        II = ABC_elti(ii);
                  end
                  
                  switch tr_typ
                    case 'P0'
                        JJ = j;
                    case 'P1'
                        JJ = ABC_eltj(jj);
                  end
                  
                  M(II,JJ) = M(II,JJ) + local_matrix(ii,jj);
                  
              end
          end
       
          % what's the convention here? hbasisx from which space?
          %local_integral = sstri_integrate(kernel,rsf_tr,rsf_ts,chi_tau,chi_t,g_tau,g_t,relation);
       
       end
    end
    
    % Creating the final matrix M

end