function [bndmesh_i,bndmesh_e] = sphSphMesh(N)
    bndmesh_i = mshSphere(N,2.5);
    bndmesh_i = bndmesh_i.translate([1 0 0]);

    bndmesh_e = mshSphere(floor(2.6*N),4);
end