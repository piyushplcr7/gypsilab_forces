
% D is the object with the floating potential / fixed charge
function [Psi,c] = solve_float_pt_ext(mesh,mesh_D,Q,order,solver,typ)

Gamma = dom(mesh,order);
Gamma_D = dom(mesh_D,order);
space = fem(mesh,typ);

switch solver
    case 'gypsi'
        % Assemble V
        Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); % 0 wave number
        V = 1/(4*pi)*integral(Gamma,Gamma,space,Gxy,space);
        V = V + 1/(4*pi)*regularize(Gamma,Gamma,space,'[1/r]',space);

    case 'ss'
        % Assemble V
        kernel = @(x,y,z) 1./vecnorm(z,2,2)/ (4*pi);
        V = panel_oriented_assembly(mesh,kernel,space,space);
end
    
    l = integral(Gamma_D,space);
    Block_matrix = [V  l;
                    l' 0];
    Psi_c = Block_matrix\[0*l; Q];
    Psi = Psi_c(1:space.ndof);
    c = Psi_c(space.ndof+1);
end