% function that returns the two traces of a vector potential 
% generated by a toroidal source centered at origin, evaluated at 
% The quadrature points specified by Gamma

function [TDA,TNA,curlAdotn,mesh_src] = getTracesTorusSource(Gamma)
    N_src = 2000;
    R0 = 2;
    r0 = .5;
    [J,mesh_src] = get_torus_source(N_src,R0,r0);
    omega_src = dom(mesh_src,3);

    [X,~] = Gamma.qud;
    
    % Computing the fields on the points X
    A = compute_vecpot(J,omega_src,X);
    curlA = compute_vecpot_curl(J,omega_src,X);
    
    % Taking traces of the fields at the points X
    normals = Gamma.qudNrm;
    TDA = A - dot(A,normals,2).*normals;
    TNA = cross(curlA,normals,2);
    curlAdotn = dot(curlA,normals,2);
end