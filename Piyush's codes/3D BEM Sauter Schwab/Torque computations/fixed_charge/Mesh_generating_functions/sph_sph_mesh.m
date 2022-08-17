function [mesh,mesh_in,mesh_out] = sph_sph_mesh(r1,r2,d,N)
    
    % distance between the centers
    dist = d+r1+r2;
    
    % Mesh for the geometry
    mesh_in = mshSphere(N,r1);
    mesh_out = mshSphere(N,r2);
    
    % Translate the outer sphere
    N_vtcs = size(mesh_out.vtx,1);
    trans = ones(N_vtcs,1) * [dist 0 0];
    mesh_out.vtx = mesh_out.vtx + trans;
    
    % Join to create the final mesh
    mesh = union(mesh_in,mesh_out);

end