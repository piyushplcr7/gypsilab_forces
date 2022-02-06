function [mesh,mesh_in,mesh_out] = cube_tor_mesh(r1,r2,L,N,seed)
    rng(seed);
    A = rand(3,3);
    [Qrot,~] = qr(A);

    % Mesh for the geometry
    mesh_out = mshTorus(N,r1,r2);
    mesh_in_vol = mshCube(N/2,L);
    mesh_in = mesh_in_vol.bnd;

    % Rotate the outer torus
    mesh_out.vtx = mesh_out.vtx * Qrot;
    
    % Translate the outer torus
    N_vtcs = size(mesh_out.vtx,1);
    trans = ones(N_vtcs,1) * [25 25 25];
    mesh_out.vtx = mesh_out.vtx + trans;
  
    % Join to create the final mesh
    mesh = union(mesh_in,mesh_out);

end