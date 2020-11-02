function [m] = mshTetra(n)

m = mshCube(1,[1 1 1]);
m = m.sub(1);
m = bnd(m);

m = m.refine(n*ones(size(m.elt,1),1));


end

