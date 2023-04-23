function [J,mesh_src] = get_torus_source(N,R0,r0)
    mesh_src = mshTorus(N,R0,r0);
    %interior_torus = @(x,y,z) sqrt( (x-x*R0./(sqrt(x.^2+y.^2))).^2 + (y-y*R0./(sqrt(x.^2+y.^2))).^2 + z.^2 )<r0;
    %torus_tangent = @(x,y,z) interior_torus(x,y,z).*[-y x z*0]./(sqrt(x.^2+y.^2));
    torus_tangent = @(x,y,z) [-y x z*0]./(sqrt(x.^2+y.^2));
    J = @(x,y,z) torus_tangent(x,y,z);
end