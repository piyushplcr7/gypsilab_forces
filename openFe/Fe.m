classdef (Abstract) Fe < fem
    % Finite Element Functions. An instance of Fe represents a vector space
    % Vh = Span{phi_1,...,phi_n} where phi_i are Finite Element functions
    % defined on a mesh. By finite element functions, we mean the following
    % conditions:
    % - There exists a set of points x_1,...,x_n on the mesh, called dofs,
    % such that phi_i(x_j) = 0 for j != i and phi_i(x_i) = 1 for all i.
    % - Each element E of the mesh is the image of the reference domain [0,1] by
    % an affine mapping chi_E. There exists B (basis) functions
    % psi_1,...,psi_B defined on S, and points X_1,...,X_B (reference dofs) of S such that
    % psi_i(X_b) = delta_{i,b}). For every f in Vh, there holds
    % $f = \sum_{E} \sum_{b = 1}^B f(x_{i(E,b)})phi_b o chi_j^{-1}$ almost
    % everywhere.
    % - Each dof x_i located on an element E coincides with chi_E(X_b) for
    % some b. If x_i is shared by several elements {e1,...,eM}, then there
    % exists b_1,...,b_m such that x_i = chi_{e_m}(X_{m}) for all m = 1..M.
    
    properties
        nb
    end
    
    methods (Abstract)
        [X,dof_jb]     =   dof(this);
        % X contains the coordinates of the dofs.
        % dof_jb(j,b) returns the index i of the dof x_i such that x_i =
        % chi_{E_j}(x_b). In other words, it returns the dof located at the
        % b-th position of the j-th element.
        [r]             =   psi_b(this,b,x);
        % input x must be a 3 containing reference coordinates.
        % the output r is a cell containing the i,j-th coordinate of psi_b(x)
        % in r{i,j}.
        [s]             =   name(this);
    end
    methods
        
        %% Access properties
        
        function[s]     =   size(this)
            r = psi_b(this,1,0);
            s = size(r);
        end
        
        function[nd]    =   ndof(this)
            [X,~] = this.dof;
            nd = size(X,1);
        end
        
        function[s]     =   dim(this)
            s = size(this.msh.elt,2);
        end
        
        function varargout = ABCD(this)
            % [A,B,C,D] = ABCD(mesh) : if mesh is a tetra mesh, A(i,:) contains the
            % coordinates of the first vtx of mesh.elts(i,:), B(i,:) the
            % second, and so on. For a particle, edge and triangle mesh, only
            % A, resp A,B, resp A,B,C, are defined.
            
            varargout = cell(1,this.dim);
            for i = 1:this.dim
                varargout{i} = this.msh.vtx(this.msh.elt(:,i),:);
            end
            for i = this.dim+1:4
                varargout{i} = NaN;
            end
        end
        
        
        %% Representation
        
        function[]      =   disp(this)
            fprintf('Finite element space of type %s on\n',name(this));
            disp(this.msh);
            fprintf('%s basis functions\n',num2str(this.nb));
            fprintf('#dofs = %s\n',num2str(size(this.elimination,2)))
        end
        
        %% Dirichlet
        
        function[P]     =   elimination(this)
            % Return elimination matrix for dirichlet 
            if isempty(this.dir)
                P = speye(this.ndof);
            else
                Wh = feval(this.name,this.dir); % Create an instance of Fe of the same class
                % On the dirichelet mesh. Wh is the vector space of functions
                % that are in Vh and not in dir(Vh).
                X = this.dof;
                [~,I] = setdiff(X,Wh.dof,'rows');
                P = sparse(I,1:length(I),1,size(X,1),length(I));
            end
            
        end
        
        %% Affine mapping to reference element:
        
        
        function[Xhat,Nabla]  =   refCoords(this,Y)
            % Evaluates the function chi_j^{-1} and its gradient at the query
            % points Y located on element Ej. More precisely, it is assumed
            % that Y is a vector of size (QxJ)x3, where J is the number of
            % element, and for each j = 1..J, Yj = Y((j-1)*Q + 1:Q,:) are points 
            % located on the element Ej. 
            % The function returns a vector Xhat and a cell Nabla such that 
            % - Xhat((j-1)*Q + 1:Q,:) are the reference coordinates of Yj
            % (thus Xhat is a (QxJ)xd where d is the dimension of the
            % reference element. 
            % - Nabla is a cell of size (d,3) and Nabal{i,j} contains the
            % derivative of the i-th component of chi_j^{-1} with respect
            % to the space variable Xj. 
            
            N = size(Y,1);
            [A,B,C,D] = ABCD(this);
            Q = N/length(this.msh);
            d = this.dim - 1; % Number of space dimensions of the reference element.
            Nabla = cell(d,3);
            switch this.dim
                case 2
                    A = repeatLines(A,Q); B = repeatLines(B,Q);
                    AY = Y-A; AB = B-A;
                    Xhat = scal3D(AY,AB)./scal3D(AB,AB);
                    NablaTmp = AB./scal3D(AB,AB);
                    Nabla{1,1} = NablaTmp(:,1); 
                    Nabla{1,2} = NablaTmp(:,2); 
                    Nabla{1,3} = NablaTmp(:,3);
                    
                case 3
                    A = repeatLines(A,Q); B = repeatLines(B,Q); C = repeatLines(C,Q);
                    [x,y,Nabla] = barycentricCoordinates(A,B,C,Y);
                    Xhat = [x,y];
                case 4
                    A = repeatLines(A,Q); B = repeatLines(B,Q); C = repeatLines(C,Q);
                    D = repeatLines(D,Q);
                    [x,y,z,Nabla] = barycentricCoordinates3(A,B,C,D,Y);
                    Xhat = [x,y,z];
            end
        end
        
        %% Linear map values at dof -> values of operator at Y.
        
        function[M]     =   uqm(this,Omega)
            Xqud = qud(Omega);
            M = this.dof2Y(Xqud);
        end
        
        function[M]     =   dof2Y(this,Y)
            % Returns the matrix M of the linear application L_Y :
            % R^{ndof} -> R^{nY} which maps the values of f \in Vh to the
            % values of op(f) at the points Y. Here, op is the operator
            % referred to in the field this.op
            % Y is a Nxd array for Fe in R^d containing points of the mesh.
            % Moreover, N is a multiple of the number of elements in the
            % mesh and is of the form [Y_1^1;...;Y_1^q;
            % Y_2^1;...;Y_2^q;...] for some q, and the points
            % Y_j^1,...,Y_j^q all lie on the element E_j.
            % If the op(phi_i) is not scalar, then M is a cell
            % and M{i} is the matrix corresponding to the i-th coordinate.
            [X,dof_jb] = this.dof;
            ndof = size(X,1);
            P = this.elimination;
            N = size(Y,1); nelt = length(this.msh); Q = N/nelt;
            assert(mod(Q,1)==0,'Y does not have the required structure');
            ind_j = cell(this.nb,1);
            ind_i = cell(this.nb,1);
            for b = 1:this.nb
                ind_j{b} = repeatLines(dof_jb(:,b),Q);
                ind_i{b} = 1:N;
            end
            switch this.opr
                case '[psi]'
                    val = feOpPsi(this,Y);
                case 'grad[psi]'
                    val = feOpGrad(this,Y);
                case 'n*[psi]'
                    val = feOpNtimes(this,Y);
                case 'x*[psi]'
                    val = feOpXtimes(this,Y);
            end
            s = size(val{1});
            if isequal(s,[1,1])
                M = sparse([],[],[],N,ndof)*P;
                for b = 1:this.nb
                    M = M + sparse(ind_i{b},ind_j{b},val{b}{1},N,ndof)*P;
                end
            else
                M = cell(s);
                for k = 1:s(1)
                    for l = 1:s(2)
                        M{k,l} = sparse([],[],[],N,ndof)*P;
                        for b = 1:this.nb
                            M{k,l} = M{k,l} ...
                                + sparse(ind_i{b},ind_j{b},val{b}{k,l},N,ndof)*P;
                        end
                    end
                end
            end
        end
        
        %% Operators. 
        % Guidelines to implement new operator: 
        % - add a case in the function dof2Y
        % - create a function function[val] = feOp<YourOperator>(this,Y) 
        % This function should return a cell val{1},...,val{nb} where nb is
        % the number of basis functions. For each b = 1..nb, val{b} is a
        % cell of size (k,l), where (k,l) is such the dimension of A[psi].
        % For example, k = 1, l = 1 for A = Id, when [psi] are scalar FE. 
        % and k = 1, l = 3 for A = grad, when [psi] are scalar FE. 
        
        function [val] = feOpPsi(this,Y)
            % operator = [psi]
            Xhat = this.refCoords(Y);
            val = cell(this.nb,1);
            for b = 1:this.nb
                val{b} = this.psi_b(b,Xhat);
            end
        end
        
        function [val] = feOpGrad(this,Y)
            % operator = grad[psi] (tangential gradient)
            try [~,d] = size(gradPsi_b(this,1,zeros(1,this.dim -1)));
                % Basis functions are P-dimensional
                % Gradient : Pxd
                % [dphi1/dx1, ..., dphi1/dxd;...; dphiP/dx1,...,dphiP/dxd])
            catch
                error(['Grad is not implemented by ' this.name])
            end
            assert(d == this.dim -1,'Error of implementation of gradient');
            [Xhat,Nabla] = this.refCoords(Y);
            val = cell(this.nb,1);
            for b = 1:this.nb
                val{b} = cellMatrixProd(this.gradPsi_b(b,Xhat),Nabla);
            end
        end
        
        function[val] = feOpNtimes(this,Y)
            Q = size(Y,1)/length(this.msh);
            nrm = this.mesh.nrm;
            nrm = repeatLines(nrm,Q);
            aux = feOpPsi(this,Y);
            val = cell(this.nb,1);
            for b = 1:this.nb
                val{b} = cell(1,3);
                val{b}{1,1} = aux{b}{1}.*nrm(:,1);
                val{b}{1,2} = aux{b}{1}.*nrm(:,2);
                val{b}{1,3} = aux{b}{1}.*nrm(:,3);
            end
        end
        
        function[val] = feOpXtimes(this,X)
            aux = feOpPsi(this,X);
            val = cell(this.nb,1);
            for b = 1:this.nb
                val{b} = cell(1,3);
                val{b}{1,1} = aux{b}{1}.*X(:,1);
                val{b}{1,2} = aux{b}{1}.*X(:,2);
                val{b}{1,3} = aux{b}{1}.*X(:,3);
            end
        end
        
        %% Operators shortcuts
        
        function[this]      =  dirichlet(this,m)
            this.dir = m;
        end
        
        function[this]      =  grad(this)
            this.opr = 'grad[psi]';
        end
        
        function[this]      = ntimes(this)
            this.opr = 'n*[psi]';
        end
        
        function[this]      = xtimes(this)
            this.opr = 'x*[psi]';
        end
        
        
        %% Value at vertices
        
        function[out]       =   valuesAtVtx(this,coords)
            % *Brilliant* idea from Gypsilab.
            Wh = P1(this.msh);
            Gamma = dom(this.msh,3);
            I = integral(Gamma,Wh,Wh);
            M = integral(Gamma,Wh,this);
            out = I \ (M*coords);
        end
    end
end

