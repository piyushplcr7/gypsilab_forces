function [curlU,Wh] = curl_NED(Vh)

if (size(Vh.msh,2)==4)
    m = Vh.msh;
    Wh = fem(m,'RWG');
    faces = m.fce;
    edges = m.edg;
    [facesSorted,~] = sort(faces.elt,2);
    faceEdges1 = [facesSorted(:,1), facesSorted(:,2)];
    faceEdges2 = [facesSorted(:,2), facesSorted(:,3)];
    faceEdges3 = [facesSorted(:,1), facesSorted(:,3)];
    [~,I1] = ismember(faceEdges1,edges.elt,'rows');
    [~,I2] = ismember(faceEdges2,edges.elt,'rows');
    [~,I3] = ismember(faceEdges3,edges.elt,'rows');
    idx = (1:m.nfce)';
    curlU = sparse(idx,I1,1,Wh.ndof,Vh.ndof) + sparse(idx,I2,1,Wh.ndof,Vh.ndof)...
        - sparse(idx,I3,1,Wh.ndof,Vh.ndof);
    
    
else
    m = Vh.msh;
    Wh = fem(m,'P0');
    gss = 1;
    Gamma = dom(m,gss);
    M = integral(Gamma,Wh,Wh);
    Fv = integral(Gamma,Wh,curl(Vh));
    curlU = spdiags(1./diag(M),0,size(M,1),size(M,2))*Fv;
    
%    % 
%     m = Vh.msh;
%     Wh = fem(m,'RWG');
%     faces = m.elt;
%     edges = m.edg;
%     [facesSorted,~] = sort(faces.elt,2);
%     faceEdges1 = [facesSorted(:,1), facesSorted(:,2)];
%     faceEdges2 = [facesSorted(:,2), facesSorted(:,3)];
%     faceEdges3 = [facesSorted(:,1), facesSorted(:,3)];
%     [~,I1] = ismember(faceEdges1,edges.elt,'rows');
%     [~,I2] = ismember(faceEdges2,edges.elt,'rows');
%     [~,I3] = ismember(faceEdges3,edges.elt,'rows');
%     idx = (1:m.nfce)';
%     curlU = sparse(idx,I1,1,Wh.ndof,Vh.ndof) + sparse(idx,I2,1,Wh.ndof,Vh.ndof)...
%         - sparse(idx,I3,1,Wh.ndof,Vh.ndof);
%     
    
end


