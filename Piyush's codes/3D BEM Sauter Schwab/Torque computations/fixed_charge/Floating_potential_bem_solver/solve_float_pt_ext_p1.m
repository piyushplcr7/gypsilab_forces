
% D is the object with the floating potential / fixed charge
function [Psi,c] = solve_float_pt_ext_p1(mesh,mesh_D,Q,order,solver)

Gamma = dom(mesh,order);
Gamma_D = dom(mesh_D,order);
S1_Gamma = fem(mesh,'P1');

switch solver
    case 'gypsi'
        % Assemble V
        Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); % 0 wave number
        V = 1/(4*pi)*integral(Gamma,Gamma,S1_Gamma,Gxy,S1_Gamma);
        V = V + 1/(4*pi)*regularize(Gamma,Gamma,S1_Gamma,'[1/r]',S1_Gamma);

    case 'ss'
        % Assemble V
        kernel = @(x,y,z) 1./vecnorm(z,2,2)/ (4*pi);
        V = panel_oriented_assembly(mesh,kernel,S1_Gamma,S1_Gamma);
end
    
    l = integral(Gamma_D,S1_Gamma);
    Block_matrix = [V  l;
                    l' 0];
    Psi_c = Block_matrix\[0*l; Q];
    Psi = Psi_c(1:S1_Gamma.ndof);
    c = Psi_c(S1_Gamma.ndof+1);
end