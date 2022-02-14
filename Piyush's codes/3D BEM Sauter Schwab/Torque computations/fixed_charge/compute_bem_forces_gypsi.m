% Function to compute BEM  forces using Gypsi implementation. 
function out = compute_bem_forces_gypsi(mesh,Psi,nu)

    % Force computation using the double layer. One arg is Psi, other is
    % Psi \nu \cdot n
    S0_Gamma = fem(mesh,'P0');
    normals = mesh.nrm;
    dofs = S0_Gamma.dof;
    Gamma = dom(mesh,3);

    GradG = cell(3,1);
    GradG{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
    GradG{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
    GradG{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);
    K = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,GradG,ntimes(S0_Gamma));
    K = K +1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'grady[1/r]',ntimes(S0_Gamma));

    % Right vectors for forces
    Psi_nu = Psi.* dot(nu(dofs),normals,2);
    out = 2*dot(Psi_nu,K*Psi);
end