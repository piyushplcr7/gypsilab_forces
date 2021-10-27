%% Getting the hard coded reference shape functions for the FEM object

% IMPORTANT: The triangular reference element is different than the one
% used in Gypsilab.

function rsfs = femRSF(fe)
    % Getting the type of FE space
    typ = fe.typ;
    mesh = fe.msh;
    
    elt_type = size(mesh.elt,2);
    
    % Segment (Is there a reference element in gypsi for line segment?)
    if elt_type == 2
        switch typ
            case 'P0'
                rsfs = cell(1,1);
                rsfs{1} = @(X) 1;
            case 'P1'
                rsfs = cell(2,1);
                rsfs{1} = @(X) 0.5 * (1+X); % support at 1
                rsfs{2} = @(X) 0.5 * (1-X); % support at -1
        end
        
    end
    % Triangle
    if elt_type == 3
       switch typ
            case 'P0'
                rsfs = cell(1,1);
                rsfs{1} = @(X) 1;
            case 'P1'
                rsfs = cell(3,1);
                %rsfs{1} = @(X) 1 - X(1) - X(2);
                %rsfs{2} = @(X) X(1);
                %rsfs{3} = @(X) X(2);
                rsfs{1} = @(X) 1 - X(1); % support at 0,0
                rsfs{2} = @(X) X(1) - X(2); % support at 1,0 
                rsfs{3} = @(X) X(2); % support at 1,1
       end
    end
end