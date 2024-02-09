function [bndmesh_i,bndmesh_e] = cubeSphMesh(N)
    bndmesh_i = getMeshCubeTranslated(N,[1 0.5 1]);

    bndmesh_e = mshSphere(40*floor(N^0.7),4);
end