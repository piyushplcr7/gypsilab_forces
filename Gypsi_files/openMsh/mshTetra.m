function [m] = mshTetra(h)

m = mshCube(1,[1 1 1]);
m = m.sub(1);
m = bnd(m);

m = m.refine(h);


end

