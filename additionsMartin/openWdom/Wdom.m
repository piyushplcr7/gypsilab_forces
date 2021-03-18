classdef Wdom < dom
    % Weighted integration domain
    
    properties
        w % int_{Wdom} f(x) w(x) dx
        dw % Gradient of w. 
        singularities = [];% Location and power of weight singularity
        Q;
        W;
        elt2qud;
        sA
        sB
        x;
        r_1x;
        L;
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
                this.singularities = sing;
            end
            this = makeQud(this);
        end
        function[this] = supplyDw(this,ddw)
            this.dw = ddw;
        end
        function[QQ,WW,eelt2qud] = qud(this)
            QQ = this.Q;
            WW = this.W;
            eelt2qud = this.elt2qud;
        end
        function[this] = makeQud(this)
            tol = 1e-8;
            g = this.gss;
            m = this.msh;
            [A,B] = m.ABCD;
            [xhat,what] = domGauss_Legendre1D(g,0,1);
            xhat = sort(xhat);
            Xhat = repmat(xhat(:),size(A,1),1);
            Arep = repeatLines(A,g);
            Brep = repeatLines(B,g);
            X0 = 0*Arep;
            X0(:,1) = Arep(:,1) + Xhat.*(Brep(:,1) - Arep(:,1));
            X0(:,2) = Arep(:,2) + Xhat.*(Brep(:,2) - Arep(:,2));
            
            omega = this.w;
            sing = this.singularities;
            singVtx = sing{1};
            singPow = sing{2};
            
            Bs = WprimitiveOnMesh(m,B,omega,singVtx,singPow,20);
            this.sB = cumsum(Bs);
            this.L = this.sB(end);
            this.sA = [0;this.sB(1:end-1)];
            Asrep = repeatLines(this.sA,g);
            Bsrep = repeatLines(this.sB,g);
            S = Xhat.*(Bsrep - Asrep);
            WW = repmat(what,size(A,1),1);
            this.W = WW.*(Bsrep(:,1) - Asrep(:,1));
            
            
            [this.Q,this.x,this.r_1x] = WfindInverses(m,omega,singVtx,singPow,S,X0,tol,this.sB);
            
            this.elt2qud = reshape(1:(g*length(this.msh)),g,length(this.msh))';
        end
        function[sY] = r_1(this,seg,Y)
            A = this.msh.vtx(this.msh.elt(seg,1),:);
            B = this.msh.vtx(this.msh.elt(seg,2),:);
            l = norm3D(B-A);
            tau = (B-A)/l;
            ssA = this.sA(seg);
            ssB = this.sB(seg);
            QQ = this.Q(this.gss*(seg-1) + 1:this.gss*seg,:);
            q = scal3D(QQ - A,tau);
            ssQ = this.r_1x(this.gss*(seg-1) + 1:this.gss*seg);
            y = scal3D(Y-A,tau);
            sY = interp1([0;q;l],[ssA;ssQ;ssB],y,'linear','extrap');
        end
        
        
    end
end

