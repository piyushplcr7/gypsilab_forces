% Function to generate mesh for cuboid 5 

function bndmesh = getMeshTetraNew(N)
    i = log(N)/log(2)-6;
    tetra_function_name = sprintf('tetra%d', i);
    tetra_function_handle = str2func(tetra_function_name);
    bndmesh = genMeshFromScript(tetra_function_handle);
    bndmesh = bndmesh.translate([2 1 3]);
end