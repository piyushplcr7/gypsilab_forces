% Helper function for Shape derivative computations

% This function computes the Neumann Trace for a Dirichlet BVP and the
% Adjoint solution for the interior Dirichlet BVP with BC given by g.
% The trial and test spaces used
% are the lowest order spaces of piecewise constant functions

function [Psi,Rho] = compute_psi_rho_int(mesh,g)
    % Definition of FEM spaces and integration rule
    S0_Gamma = fem(mesh,'P0');
    order = 3;
    Gamma = dom(mesh,order);
    
    % Getting the Single Layer matrix V using Gypsilab implementation
    Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); % 0 wave number
    V = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,Gxy,S0_Gamma);
    V = V + 1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'[1/r]',S0_Gamma);
    
    % Getting the Double Layer matrix K using Gypsilab implementation
    GradG = cell(3,1);
    GradG{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
    GradG{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
    GradG{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);
    K = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,GradG,ntimes(S0_Gamma));
    K = K +1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'grady[1/r]',ntimes(S0_Gamma));
    
    % Defining the mass matrix M
    M = integral(Gamma,S0_Gamma,S0_Gamma);
    
    % Constructing \int_{Gamma} g phi dS for all test functions phi
    g_N = integral(Gamma,S0_Gamma,g);
    
    % Solving exterior problem, M\g_N = coeffs of g in terms of S0_Gamma
    Psi = V\((0.5 * M + K)* (M\g_N));
    
    % Solving the adjoint problem to get the adjoint solution
    Rho = V\(-0.5 * g_N);
end