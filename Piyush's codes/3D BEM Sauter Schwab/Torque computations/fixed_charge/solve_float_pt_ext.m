
% D is the object with the floating potential / fixed charge
function [Psi,c] = solve_float_pt_ext(mesh,mesh_D,Q,order,solver)

Gamma = dom(mesh,order);
Gamma_D = dom(mesh_D,order);
S0_Gamma = fem(mesh,'P0');

switch solver
    case 'gypsi'
        % Assemble V
        Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); % 0 wave number
        V = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,Gxy,S0_Gamma);
        V = V + 1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'[1/r]',S0_Gamma);

    case 'ss'
        % Assemble V
        kernel = @(x,y,z) 1./vecnorm(z,2,2)/ (4*pi);
        V = panel_oriented_assembly(mesh,kernel,S0_Gamma,S0_Gamma);
end
    
    l = integral(Gamma_D,S0_Gamma);
    Block_matrix = [V  l;
                    l' 0];
    Psi_c = Block_matrix\[0*l; Q];
    Psi = Psi_c(1:S0_Gamma.ndof);
    c = Psi_c(S0_Gamma.ndof+1);
end