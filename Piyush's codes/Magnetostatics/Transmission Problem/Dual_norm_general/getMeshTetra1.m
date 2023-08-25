% Function to generate mesh for tetra

function bndmesh = getMeshTetra1(N)
    % N from 128 to 2048  
    i = log(N)/log(2);
    % For lowest N = 128 (2^7) i should be 1
    i = i-6;
    vals = [1 1/2 1/4 1/7.9 1/15.9 1/31.9];
    bndmesh = meshSymTetra;
    bndmesh = bndmesh.translate([2 1 3]);
    bndmesh = bndmesh.refine(vals(i));
end