POS = msh.POS;
TETS = msh.TETS;
TRGLS = msh.TRIANGLES;
clear msh;
trmesh = msh(POS,TRGLS(:,1:3));
tetmesh = msh(POS,TETS(:,1:4));