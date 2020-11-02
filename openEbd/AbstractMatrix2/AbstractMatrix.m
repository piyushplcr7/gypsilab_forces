classdef AbstractMatrix
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        concretePart;
        abstractPart;
        N1;
        N2;
    end
    
    methods
        function[this] = AbstractMatrix(cP,aP,NN1,NN2)
            if nargin == 0
                this.concretePart = 0;
                this.abstractPart = @(x)(0);
                this.N1 = 1;
                this.N2 = 1;
            elseif nargin == 1
                if isa(cP,'AbstractMatrix')
                    this.concretePart = cP.concretePart;
                    this.abstractPart = cP.abstractPart;
                    this.N1 = size(cP,1);
                    this.N2 = size(cP,2);
                elseif(isa(cP,'double'))
                    this.concretePart = cP;
                    this.abstractPart = [];
                    this.N1 = size(cP,1);
                    this.N2 = size(cP,2);
                end
            elseif nargin ==2
                assert(~isempty(cP),'cannot guess matrix size, specify N1 and N2');
                this.concretePart = cP;
                this.abstractPart = aP;
                this.N1 = size(cP,1);
                this.N2 = size(cP,2);
            elseif nargin == 4
                if isempty(cP)
                    cP = sparse(NN1,NN2);
                end
                this.concretePart = cP;
                this.abstractPart = aP;
                this.N1 = NN1;
                this.N2 = NN2;
                assert(isequal([NN1 NN2],size(cP)));
            end
        end
        function[out] = size(this,dim)
            if nargin==1
                out = [this.N1,this.N2];
            else
                if dim ==1
                    out = this.N1;
                else
                    out = this.N2;
                end
            end
        end
        function[b] = isConcrete(this)
            b = isempty(this.abstractPart);
        end
        function[mv] = abstractMV(this)
            if this.isConcrete
                mv = @(x)(zeros(this.N1,1));
            else
                mv = this.abstractPart;
            end
        end
        function[this_ij] = elem(this,i,j)
            zj = zeros(size(this,2),1);
            zi = zeros(size(this,1),1);
            ei = zi; ei(i) = 1;
            ej = zj; ej(j) = 1;
            this_ij = ei'*(this*ej);
        end
        function[C] = plus(A,B)
            if ~isa(A,'AbstractMatrix')
                C = plus(B,A);
                return
            end
            if and(isa(B,'double'),isscalar(B))
                C = A+B*AbstractMatrix.ones(size(A));
            elseif and(isa(A,'double'),isscalar(A))
                C = plus(B,A);
            else
                assert(isequal(size(A),size(B)),'Error using  + (Abstract) Matrix dimensions must agree.');
                C = A;
                A = AbstractMatrix(A);
                B = AbstractMatrix(B);                
                C.concretePart = A.concretePart + B.concretePart;
                if ~and(A.isConcrete,B.isConcrete)
                    mvA = A.abstractMV;
                    mvB = B.abstractMV;
                    C.abstractPart = @(x)(mvA(x) + mvB(x));
                else
                    C.abstractPart = [];
                end
            end
        end
        function[C] = uminus(A)
            C = A;
            if isa(A,'AbstractMatrix')
                C.concretePart = -A.concretePart;
                if ~A.isConcrete
                    mvA = A.abstractPart;
                    C.abstractPart = @(x)(-mvA(x));
                end
            else
                error('this statement should not be reachable');
            end
        end
        function[C] = minus(A,B)
            C = plus(A,-B);
        end
        function[C] = times(A,B)
            if or(isequal(size(A),[1,1]),isequal(size(B),[1,1]))
                C = A*B;
            else
                assert(isequal(size(A),size(B)),'Matrix dimensions do not agree');
                yn = input('You are about to convert matrix to concrete matrix, proceed (y/n) ? \n');
                if isequal(yn,'y')
                    C = concrete(A).*concrete(B);
                else
                    error('Computation interrupted')
                end
            end
        end
        function[C] = mtimes(A,B)
            assert(or(isa(B,'double'),isa(B,'AbstractMatrix')),sprintf(...
                '2nd arg is of type %s but should be of type ''AbstractMatrix'' or ''double''',...
                class(B)));
            assert(or(isa(A,'double'),isa(A,'AbstractMatrix')),sprintf(...
                '1st arg is of type %s but should be of type ''AbstractMatrix'' or ''double''',...
                class(B)));
            assert(or(or(size(A,2)==size(B,1),isscalar(A)),isscalar(B)));
            NN1 = size(A,1);
            NN2 = size(B,2);
            if or(isscalar(A),isscalar(B))
                if isscalar(A)
                    assert(isa(B,'AbstractMatrix'));
                    C = B;
                    C.concretePart = C.concretePart*A;
                    if ~C.isConcrete
                        mvC = C.abstractPart;
                        C.abstractPart = @(x)(A*mvC(x));
                    end
                else
                    % B is scalar
                    C = mtimes(B,A);
                end
                return
            end
            if isa(A,'double')
                A = AbstractMatrix(A);
            else
                assert(isa(A,'AbstractMatrix'));                
            end
            cA = A.concretePart;
            mvA = A.abstractMV;
            if and(isa(B,'double'),size(B,2)==1)                
                C = cA*B + mvA(B);
                return
            elseif isa(B,'double')
                B = AbstractMatrix(B);
            else
                assert(isa(B,'AbstractMatrix'));    
            end
            cB = B.concretePart;
            cC = cA*cB;
            if B.isConcrete
                mvC = @(x)(mvA(cB*x));
            elseif A.isConcrete
                mvB = B.abstractPart;
                mvC = @(x)(cA*mvB(x));
            else
                mvB = B.abstractPart;
                mvTemp = @(x)(mvA(x) + cA*x);
                mvC = @(x)(mvTemp(mvB(x) + cB*x) - cA*cB*x);
            end
            if and(A.isConcrete,B.isConcrete)
                C = AbstractMatrix(cC,[],NN1,NN2);
            else                
                C = AbstractMatrix(cC,mvC,NN1,NN2);
            end
        end
        function[A] = concrete(A)
            if A.isConcrete
                % nothing to do;
            else
                mv = A.abstractMV;
                I = eye(A.N2);
                for i = 1:A.N2
                    ei = I(:,i);
                    A.concretePart(:,i) = A.concretePart(:,i) + mv(ei);
                end
                A.abstractPart = [];
            end
        end
        function[A] = cP(A)
            A.abstractPart = []; 
        end
        function[A] = aP(A)
            A.concretePart = [];
        end
        function[C] = full(A)
            A = A.concrete;
            B = A.concretePart;
            C = full(B);
        end
        function[this] = real(this)
            this.concretePart = real(this.concretePart);
            if ~this.isConcrete
                mv = this.abstractPart;
                this.abstractPart = @(x)(real(mv(x)));
            end
        end
        function[x,flag,relres,iter,resvec] = mldivide(this,b,restart,tol,maxit,M1,M2,x0)  
            if all([nargin==2,nargout==1,this.isConcrete])
                A = this.concretePart;
                x = A\b;
            else
                if ~exist('restart','var')
                    restart = [];
                end
                if ~exist('tol','var')
                    tol = [];
                end
                if ~exist('maxit','var')
                    maxit = [];
                end
                if ~exist('M1','var')
                    M1 = [];
                end
                if ~exist('M2','var')
                    M2 = [];
                end
                if ~exist('x0','var')
                    x0 = [];
                end
                [x,flag,relres,iter,resvec] = gmres(@(x)(this*x),b,restart,tol,maxit,M1,M2,x0);                
            end
        end
        function[x,flag,relres,iter,resvec] = gmres(this,b,varargin)
            [x,flag,relres,iter,resvec] = mldivide(this,b,varargin{:});
        end
        function[] = spy(this)
            spy(this.concretePart);
        end
        function[out] = nnz(this)
            out = nnz(this.concretePart);
        end
        function[str] = concreteMem(this)
            bytes = max(nnz(this),1) * (12) + (this.N1+1)*4;
            str = Bytes2str(bytes);
        end
    end
    
    methods (Static)
        function[out] = ones(Nx,Ny)
            if nargin == 1
                if length(Nx)==1
                    Ny = Nx;
                else
                    Ny = Nx(2);
                    Nx = Nx(1);
                end
            end
            cP = sparse([],[],[],Nx,Ny);
            aP = @(x)(ones(Nx,1)*sum(x));
            out = AbstractMatrix(cP,aP,Nx,Ny);
        end
        function[out] = rank1(u,v)
            % Returns the matrix given by uv'
            u = u(:);
            v = v(:);
            
            aP = @(x)((sum(v.*x))*u);
            N1 = length(u);
            N2 = length(v);
            cP = sparse([],[],[],N1,N2);
            out = AbstractMatrix(cP,aP,N1,N2);
        end
        function[out] = spdiag(vec)
            vec = vec(:);
            N = length(vec);
            out = AbstractMatrix(spdiags(vec,0,N,N));
        end
        function[out] = eye(m,storageOpt)
            if nargin == 1 || isequal(storageOpt,'abstract')
                cP = sparse([],[],[],m,m);
                aP = @(x)(x);
            else
                if isequal(storageOpt,'concrete')
                    cP = speye(m);
                    aP = [];
                else
                    error('2nd argument not recognized')
                end
            end
            out = AbstractMatrix(cP,aP,m,m);
        end
    end
    
end

