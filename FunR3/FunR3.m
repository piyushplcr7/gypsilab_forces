classdef FunR3 < handle
    % Function R^2 -> R
    properties
        f
    end
    
    methods
        
        %% Class constructor
        
        
        function[this] = FunR3(arg)
            if nargin == 0
                this.f = @(Z)(sparse(size(Z,1),1)); % Identical 0 func
            else
                if isa(arg,'function_handle')
                    assert(isequal(size(arg(zeros(3,3))),[3,1]));
                    this.f = arg;
                elseif isa(arg,'double')
                    if isempty(arg)
                        this.f = @(Z)(sparse(size(Z,1),1));
                    else
                        assert(isscalar(arg));
                        this.f = @(Z)(0*Z(:,1) + arg);
                    end
                end
            end
        end
        
        %% Display
        
        function[] = disp(~)
            fprintf('R3toCfunc object \n');
        end
        
        
        %% calling method
        
        
        function[v] = subsref(this,Z)
            assert(length(Z.subs)==1,'too many input arguments');
            assert(size(Z.subs{1},2)==3,...
                'function must be evaluated on Nx3 args');
            v = this.f(Z.subs{1});
        end
        
        %% Algrebraic operations 
        
        function[C] = plus(A,B)
            if isa(A,'double')
                C = FunR3(@(Z)(A + B.f(Z)));
            elseif isa(B,'double')
                C = FunR3(@(Z)(A.f(Z) + B));
            else
                C = FunR3(@(Z)(A.f(Z) + B.f(Z)));
            end
        end
        function[C] = minus(A,B)
            C = plus(A,-1*B);
        end
        function[C] = uminus(A)
            C = -1*A;
        end
        function[C] = mtimes(A,B)
            if and(isa(A,'double'),isscalar(A))
                C = FunR3(@(Z)(A*B.f(Z)));
            elseif and(isa(B,'double'),isscalar(B))
                C = FunR3(@(Z)(A.f(Z)*B));
            else
                assert(isequal(class(A),'FunR3'));
                assert(isequal(class(B),'FunR3'));
                C = FunR3(@(Z)(A.f(Z).*B.f(Z)));
            end
        end
        function[C] = mpower(A,B)
            if and(isa(A,'double'),isscalar(A))
                C = FunR3(@(Z)(A.^(B.f(Z))));
                
            elseif and(isa(B,'double'),isscalar(B))
                C = FunR3(@(Z)(A.f(Z).^B));
            else
                assert(isequal(class(A),'FunR3'));
                assert(isequal(class(B),'FunR3'));
            end
        end
        function[C] = times(A,B)
            C = mtimes(A,B);
            
        end
        function[C] = rdivide(A,B)
            if isa(A,'double')
                C = FunR3(@(Z)(A./B.f(Z)));     
            elseif isa(B,'double')
                C = FunR3(@(Z)(A.f(Z)./B));
            else
                C = FunR3(@(Z)(A.f(Z)./B.f(Z)));
            end
        end
        function[C] = mrdivide(A,B)
            C = rdivide(A,B);
        end
        function[B] = applyFun(A,arbitraryFun)
            B = FunR3(@(Z)(arbitraryFun(A.f(Z))));
        end
        function[B] = conj(A)
            B = applyFun(A,@conj);
        end
        function[B] = exp(A)
            B = applyFun(A,@exp);
        end
        function[B] = sin(A)
            B = applyFun(A,@sin);
        end
        function[B] = cos(A)
            B = applyFun(A,@cos);
        end
        function[B] = sqrt(A)
            B = applyFun(A,@sqrt);
        end
        function[B] = ln(A)
            B = applyFun(A,@log);
        end
        function[B] = log(A)
            B = applyFun(A,@log);
        end
        function[C] = power(A,n)
            C = mpower(A,n);
        end
        function[C] = abs(A)
            C = applyFun(A,@abs);
        end
        function[C] = atan2(Y,X)
            C = FunR3(@(Z)(atan2(Y.f(Z),X.f(Z))));
        end
        function[B] = asin(X)
            B = applyFun(X,@asin);
        end
        function[B] = acos(X)
            B = applyFun(X,@acos);
        end
        function[B] = dilog(A)
            B = applyFun(A,@dilog);
        end
        
        
    end
    methods (Static)
        function[B] = X
            B = FunR3(@(Z)(Z(:,1)));
        end
        function[B] = Y
            B = FunR3(@(Z)(Z(:,2)));
        end
        function[B] = Z
            B = FunR3(@(Z)(Z(:,3)));
        end
        function[A,B,C] = XYZ
            A = FunR3.X; B = FunR3.Y; C = FunR3.Z;
        end
    end
end

