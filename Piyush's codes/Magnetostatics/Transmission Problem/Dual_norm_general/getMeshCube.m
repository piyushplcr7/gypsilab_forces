% Function to generate mesh for cube
% L = 2*[1 1 1];
%     T = [5 5 3];

function bndmesh = getMeshCube(N)
    L = 2*[1 1 1];
    T = [5 5 3];
    
    mesh = mshCube(N,L);
    mesh = mesh.translate(T);
    bndmesh = mesh.bnd;
end