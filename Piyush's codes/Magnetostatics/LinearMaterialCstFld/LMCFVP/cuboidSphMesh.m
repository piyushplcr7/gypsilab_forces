function [bndmesh_i,bndmesh_e] = cuboidSphMesh(N)
    bndmesh_i = getMeshCuboid5Translated(floor(N/3),[1 0.5 1]);

    bndmesh_e = mshSphere(40*floor(N^0.7),4);
end