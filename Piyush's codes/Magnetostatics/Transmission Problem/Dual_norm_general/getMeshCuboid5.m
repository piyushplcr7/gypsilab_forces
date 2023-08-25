% Function to generate mesh for cuboid 5 

function bndmesh = getMeshCuboid5(N)
    L = [3 1 1];
    T = [2 1 3];
    
    mesh = mshCube(N,L);
    mesh = mesh.translate(T);
    bndmesh = mesh.bnd;
end