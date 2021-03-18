classdef P1 < Fe
    
    methods
        
        %% Constructor
        
        function[this] = P1(m)
            this.typ = 'P1';
            this.opr = '[psi]';
            this.msh = m;
            this.nb = this.dim;
        end
        
        %% Must implement:
        
        function[s]    = name(~)
            s = 'P1';
        end
        
        function[X,elt2dof] = dof(this)
            X = this.msh.vtx;
            elt2dof = this.msh.elt;
        end
        
        function[out] = psi_b(this,b,X)
            switch this.dim
                case 2 % Edge mesh
                    x = X; % X is one-dimensional
                    switch b
                        case 1
                            r = 1 - x;
                        case 2
                            r = x;
                    end
                case 3 % Triangle mesh
                    x = X(:,1); y = X(:,2);
                    switch b
                        case 1
                            r = 1 - x - y;
                        case 2
                            r = x;
                        case 3
                            r = y;
                    end
                case 4 % Tetrahedral mesh
                    x = X(:,1); y = X(:,2); z = X(:,3);
                    switch b
                        case 1
                            r = 1 - x - y - z;
                        case 2
                            r = x;
                        case 3
                            r = y;
                        case 4
                            r = z;
                    end
            end
            out = {r};
        end
        
        
        %% Gradient
        
        function[r] = gradPsi_b(this,b,X)
            switch this.dim
                case 2 % Edge mesh
                    x = X;
                    switch b
                        case 1
                            r = {0*x - 1};
                        case 2
                            r = {0*x + 1};
                    end
                case 3 % Triangular mesh
                    r = cell(1,2);
                    x = X(:,1); y = X(:,2);
                    switch b
                        case 1
                            r{1,1} = 0*x - 1; 
                            r{1,2} = 0*y - 1; 
                        case 2
                            r{1,1} = 0*x + 1;
                            r{1,2} = 0*y;
                        case 3
                            r{1,1} = 0*x;
                            r{1,2} = 0*y + 1;
                    end
                    
            end
        end
    end
    
end