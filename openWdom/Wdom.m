classdef Wdom < dom
    % Weighted integration domain
    
    properties
        w % int_{Wdom} f(x) w(x) dx
        singularities = [];% Location and power of weight singularity
    end
    
    methods
        %% Constructor
        function[this] = Wdom(m,g,w,sing)
            assert(size(m.elt,2)==2,'Wdom only supported for edge meshes');
            this.msh = m; this.gss = g; 
            if nargin >= 3
                this.w = w;
            end
            if nargin == 4
                if iscell(sing)
                    singVtx = sing{1};
                    [~,I] = intersect(m.vtx,singVtx,'rows');
                    sing = [I(:),sing{2}(:)];
                end
                this.singularities = sing;
            end
        end
        function[Q,W,elt2qud] = qud(this)
            g = this.gss;
            m = this.msh;
            nelt = length(m);
            vtx = m.vtx; elt = m.elt;
            A = vtx(elt(:,1),:); B = vtx(elt(:,2),:);
            [xhat,what] = domGauss_Legendre1D(g,0,1);
            xhat = sort(xhat);
            Xhat = repmat(xhat(:),size(A,1),1);
            What = repmat(what(:),size(A,1),1);
            A = repeatLines(A,g);
            B = repeatLines(B,g);
            Q = 0*A;
            Q(:,1) = A(:,1) + Xhat.*(B(:,1) - A(:,1));
            Q(:,2) = A(:,2) + Xhat.*(B(:,2) - A(:,2));
            W = What.*repeatLines(m.ndv,g);
            for i = 1:size(this.singularities,1)
                sing_i = this.singularities(i,:);
                vtx_i = sing_i(1); alpha_i = sing_i(2);
                beta_i = 1./(1 + alpha_i);
                eltSort = sort(elt,2);
                
                elts = find(ismember(eltSort,vtx_i));
                for j = 1:length(elts)
                    eltj = elts(j,:);
                    indj = g*(eltj-1) + (1:g);
                    if eltj(1) == vtx_i
                        uhat = xhat;
                        u = xhat.^beta_i;
                    else
                        uhat = 1 - xhat;
                        u = 1 - (1 - xhat).^beta_i;
                    end
                    
                    Q(indj,1) = A(indj,1) + u.*(B(indj,1) - A(indj,1));
                    Q(indj,2) = A(indj,2) + u.*(B(indj,2) - A(indj,2));
                    
                    W(indj) = W(indj).*beta_i.*uhat.^(beta_i - 1);
                end
            end
            W = W.*this.w(Q);
            elt2qud = reshape(1:(g*nelt),g,nelt)';
        end
        
    end
end

