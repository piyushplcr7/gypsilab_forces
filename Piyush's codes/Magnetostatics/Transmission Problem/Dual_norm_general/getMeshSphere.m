% Function to generate mesh for cube
% R = 1
%     T = [5 5 3];

function bndmesh = getMeshSphere(N)
    T = [5 5 3];
    bndmesh = mshSphere(N,1);
    bndmesh = bndmesh.translate(T);
end