% Function to compute BEM  forces using Gypsi implementation of Double Layer Galerkin matrix. 
function out = compute_bem_forces_gypsi(mesh,Psi,nu)

    % Force computation using the double layer. One arg is Psi, other is
    % Psi \nu \cdot n
    S0_Gamma = fem(mesh,'P0');
    normals = mesh.nrm;
    qudnrms = repelem(normals,3,1);
    dofs = S0_Gamma.dof;
    Gamma = dom(mesh,3);
    [X,W,elt2qud] = Gamma.qud;
    uqmat = S0_Gamma.uqm(Gamma);
    nudotn = dot(nu(X),qudnrms,2);

    GradG = cell(3,1);
    GradG{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
    GradG{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
    GradG{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);
    K = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,GradG,ntimes(S0_Gamma));
    K = K +1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'grady[1/r]',ntimes(S0_Gamma));

    % Right vectors for forces

    % Nodal interpolation kinda thing
    %Psi_nu = Psi.* dot(nu(dofs),normals,2);

    % L2 projection
    % Mass matrix
    M = integral(Gamma,S0_Gamma,S0_Gamma);

    % Rhs vector
    % Constructing the relevant matrix first
    N = S0_Gamma.ndof;
    rhsmat = spalloc(N,N,3*N);
 
    for i=1:N
        for j=1:N
            uqi = uqmat(:,i);
            uqj = uqmat(:,j);
            if normest(uqi-uqj)==0
                rhsmat(i,j) = sum(W.*uqi.*uqj.*nudotn);
            end
        end
    end
    Psi_nu = M\(rhsmat * Psi);
    out = dot(Psi_nu,K*Psi);
end