% This function computes the force and torque on object D
function [torque,force] = compute_mst_forces(mesh,Xcg,Psi)
    % Creating the FEM spaces
    S0_Gamma = fem(mesh,'P0');

    % Creating the r X n vector
    dofs = S0_Gamma.dof;
    Xvec = dofs-ones(size(dofs,1),1)*Xcg;
    normals = mesh.nrm;
    torque = 0.5 * sum( mesh.ndv .* Psi.^2 .* cross(Xvec,normals) ,1);
    force = 0.5 * sum( mesh.ndv .* Psi.^2 .* normals ,1);
end