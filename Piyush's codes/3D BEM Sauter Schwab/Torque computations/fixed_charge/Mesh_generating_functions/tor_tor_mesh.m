function [mesh,mesh_in,mesh_out] = tor_tor_mesh(r1,r2,N,seed)
    rng(seed);
    A = rand(3,3);
    [Qrot,~] = qr(A);
    
    % distance between the centers
    dist = 3*(r1+r2);
    
    % Mesh for the geometry
    mesh_in = mshTorus(N,r1,r2);
    mesh_out = mesh_in;
    
    % Translate the outer torus
    N_vtcs = size(mesh_out.vtx,1);
    trans = ones(N_vtcs,1) * [dist 0 0];
    mesh_out.vtx = mesh_out.vtx + trans;
    
    % Rotate the inner torus
    mesh_in.vtx = mesh_in.vtx * Qrot;
    
    % Join to create the final mesh
    mesh = union(mesh_in,mesh_out);

end