function [bndmesh_i,bndmesh_e] = tetraSphMesh(N)
    bndmesh_i = getMeshTetraNew(N);

    bndmesh_e = mshSphere(floor(N^2/287),2);
        bndmesh_e = bndmesh_e.translate([2 1 3]);
        bndmesh_e = bndmesh_e.translate([0.3 0.5 0.1]);
end