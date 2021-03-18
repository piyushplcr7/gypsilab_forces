function [piTU,Wh] = tangentialTrace_NED(Vh)



m = Vh.msh;
dm = m.bnd;
Wh = fem(dm,'NED');
ndofV = Vh.ndof;
ndofW = Wh.ndof;
mshedg = m.edg;
dmshedg = dm.edg;
[~,I] = ismembertol(dm.vtx,m.vtx,'ByRows',true);
idx = (1:ndofW)';
[~,jdx] = ismember(I(dmshedg.elt),mshedg.elt,'rows');
piTU = -sparse(idx,jdx,1,ndofW,ndofV);


end

