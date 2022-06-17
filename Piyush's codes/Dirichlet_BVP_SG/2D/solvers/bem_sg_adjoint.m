% This function is used to solve for the Adjoint solution 
% inputs : boundary mesh, function for boundary data g 

function p = bem_sg_adjoint(mesh,g)
% Setting up the discrete system

% Quadrature order
order = 3;
% Creating the integration domain
Gamma = dom(mesh, order);

% Defining the BEM spaces
S1_Gamma = fem(mesh,'P1'); % Space for discretizing g
S0_Gamma = fem(mesh,'P0'); % Usual space for Neumann trace space?

% Defining the kernel for SL
Gxy = @(X,Y)femGreenKernel(X,Y,'[log(r)]',0); % 0 wave number
% Defining the kernel for DL
GradGn = cell(3,1);
GradGn{1} = @(X,Y)femGreenKernel(X,Y,'grady[log(r)]1',0);
GradGn{2} = @(X,Y)femGreenKernel(X,Y,'grady[log(r)]2',0);
GradGn{3} = @(X,Y)femGreenKernel(X,Y,'grady[log(r)]3',0);

% LHS
% Evaluating the bilinear form for the Single Layer BIO
V = -1/(2*pi)*integral(Gamma,Gamma,S0_Gamma,Gxy,S0_Gamma);
V = V + -1/(2*pi)*regularize(Gamma,Gamma,S0_Gamma,'[log(r)]',S0_Gamma);

% RHS
l = -0.5 * integral(Gamma,S0_Gamma,S1_Gamma);

% A clean way to get the coefficients g_N such that when combined with the
% basis elements in trial_space_g, we get an approximation of g?
g_N = g(S1_Gamma.dof);

% The linear system is V p = l * g_N;

% Neumann trace solution:
p = V \ (l * g_N );
    
end