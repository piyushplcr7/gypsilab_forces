classdef P0 < Fe

    methods
        function[this] = P0(m)
            this.mesh = m;
            this.nb = 1;
            this.ndof = m.nelt;
        end
        function[s]    = name(~)
            s = 'P0';
        end
        function[X,elt2dof] = dof(this)
            X = this.mesh.mid;
            elt2dof = (1:this.mesh.nelt)';
        end
        function[r] = psi_b(~,~,x)
            r{1} = 0*x + 1;
        end
    end
    
end