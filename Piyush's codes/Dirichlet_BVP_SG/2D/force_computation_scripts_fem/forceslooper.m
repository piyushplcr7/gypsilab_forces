fileID = fopen('results.txt','w');
fprintf(fileID,'#Fx_bnd, Fy_bnd; Fx_vol, Fy_vol \n');
fclose(fileID);

sqkite_M_20;
vtx = msh.POS;
elt = msh.TRIANGLES(:,1:3);
clear msh;
mesh = msh(vtx,elt);
forces(mesh);
clear;

sqkite_M_30;
vtx = msh.POS;
elt = msh.TRIANGLES(:,1:3);
clear msh;
mesh = msh(vtx,elt);
forces(mesh);
clear;

sqkite_M_40;
vtx = msh.POS;
elt = msh.TRIANGLES(:,1:3);
clear msh;
mesh = msh(vtx,elt);
forces(mesh);
clear;

sqkite_M_50;
vtx = msh.POS;
elt = msh.TRIANGLES(:,1:3);
clear msh;
mesh = msh(vtx,elt);
forces(mesh);
clear;

sqkite_M_60;
vtx = msh.POS;
elt = msh.TRIANGLES(:,1:3);
clear msh;
mesh = msh(vtx,elt);
forces(mesh);
clear;

sqkite_M_70;
vtx = msh.POS;
elt = msh.TRIANGLES(:,1:3);
clear msh;
mesh = msh(vtx,elt);
forces(mesh);
clear;

sqkite_M_80;
vtx = msh.POS;
elt = msh.TRIANGLES(:,1:3);
clear msh;
mesh = msh(vtx,elt);
forces(mesh);
clear;

sqkite_M_90;
vtx = msh.POS;
elt = msh.TRIANGLES(:,1:3);
clear msh;
mesh = msh(vtx,elt);
forces(mesh);
clear;