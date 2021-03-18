function[pinU,Wh] = normalTrace_RT(Vh)


tol = 1e-12;

[Xdof,elt2dof] = Vh.dof;
m = Vh.msh;
dm = m.bnd;

Wh = fem(dm,'P0');

Ndof = size(Xdof,1);
normals = dm.nrm;

fce = m.fce;
faceCenters = unique(fce.ctr,'rows');
Ngss = size(Vh.msh.elt,2);
Ife     = (1:size(Vh.msh.elt,1))';

Ndv = m.ndv;
Nbas = size(elt2dof,2);
idx = zeros(m.nelt,Nbas,Ngss);
jdx = zeros(m.nelt,Nbas,Ngss);
for i = 1:Nbas
    for j = 1:Ngss
        idx(:,i,j) = j + Nbas*(0:m.nelt-1)';
        jdx(:,i,j) = elt2dof(Ife,i);
    end
end


if size(m.elt,2) == 3
    % Triangular mesh
    nface = 3*m.nelt;
    M = cell(1,3);
    for i = 1:Nbas
        % Indices of neighbors
        ip1 = mod(i,3)+1;
        ip2 = mod(ip1,3)+1;
        faceCenters(i:Nbas:nface,:) = 1/2*...
            (m.vtx(m.elt(:,ip1),:)...
            + m.vtx(m.elt(:,ip2),:));
    end
    for n = 1:3
        % Initialization
        val = zeros(m.nelt,Nbas,Ngss);
        
        % For each basis function
        for i = 1:Nbas
            % Indices of neighbors
            ip1 = mod(i,3)+1;
            ip2 = mod(ip1,3)+1;
            
            % Flux through the edge
            flux = (2*(m.elt(:,ip1) < m.elt(:,ip2))-1)./(2*Ndv);
            
            % For each integration point
            for j = 1:Ngss
                val(:,i,j) = flux.*(faceCenters(j:Ngss:nface,n) - m.vtx(m.elt(:,i),n));
            end
        end
        
        % Operator
        M{n} = sparse(idx(:),jdx(:),val(:),3*m.nelt,Ndof);
    end
    dofsP0 = Wh.dof;
    [~,I] = ismembertol(dofsP0,faceCenters,tol,'ByRows',true);
    idx = (1:size(dofsP0,1))';
    jdx = I;
    P = sparse(idx,jdx,1,size(dofsP0,1),nface);
    pinU = zeros(size(dofsP0,1),1);
    for n =1:3
        Nn = spdiags(normals(:,n),0,Wh.ndof,Wh.ndof);
        pinU = pinU + Nn*(P*M{n});
    end
    
else
    dofsP0 = Wh.dof;
    idx = (1:size(dofsP0,1))';
    [~,jdx] = ismembertol(dofsP0,fce.ctr,'ByRows',true);
    P = sparse(idx,jdx,1,Wh.ndof,Vh.ndof);
    normals = dm.nrm;
    facesDm = dm.fce;
    fceSorted = sort(facesDm.elt,2);
    naturalOrientation = cross(facesDm.vtx(fceSorted(:,2),:) - facesDm.vtx(fceSorted(:,1),:),...
        facesDm.vtx(fceSorted(:,3),:) - facesDm.vtx(fceSorted(:,1),:),2);
    cst = -1./dm.ndv.*sign(dot(normals,naturalOrientation,2));
    pinU = sparse(idx,idx,cst,Wh.ndof,Wh.ndof)*P;
    
end



end