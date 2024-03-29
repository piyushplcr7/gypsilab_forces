function [dof,elt2dof,restrictions] = dofMS(Vh)

typ = Vh.typ(7:end);
ms = Vh.msh;
pans = ms.panels;
restrictions = cell(ms.npanels,1);
elt2dof= cell(ms.npanels,1);
X = cell(ms.npanels,1);
if strcmp(typ,'P1')
    dof = [];
    for i = 1:ms.npanels
        femi = fem(pans{i},typ);
        [X{i},elt2dof{i}] = femi.dof;
        dof = [dof; [X{i},i*ones(size(X{i},1),1)]];
    end
    bound = ms.bnd;
    [~,jdx] = ismembertol(dof(:,1:3),bound.vtx,'ByRows',true);
    dof = dof(~jdx,:);
    dofBound = zeros(bound.nvtx,4);
    dofBound(:,1:3) = bound.vtx; dofBound(:,4) = 0;
    dof = [dofBound;dof];
    for i = 1:ms.npanels
        Xi = [X{i},ones(size(X{i},1),1)*i];
        [~,jdx] = ismembertol(Xi,dof,'ByRows',true);
        idx1 = find(jdx~=0);
        jdx1 = jdx(jdx~=0);
        Xi(:,4) = 0;
        [~,jdx] = ismembertol(Xi,dof,'ByRows',true);
        idx2 = find(jdx~=0);
        jdx2 = jdx(jdx~=0);
        idx = [idx1;idx2];
        jdx = [jdx1;jdx2];
        elt2dof{i}(idx,:) = jdx(elt2dof{i}(idx,:));
        restrictions{i} = sparse(idx,jdx,1,size(X{i},1),size(dof,1));
    end
    
else
    error('unavailable case')
end






end

