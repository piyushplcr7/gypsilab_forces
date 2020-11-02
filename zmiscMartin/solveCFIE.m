function [lambda,niter,resvec,tSolve] = solveCFIE(m,k,u,varargin)

p = inputParser;
p.addOptional('niter',size(m.vtx,1));
p.parse(varargin{:});

Gamma = dom(m,3);
Vh = fem(m,'P1');
[X,Y,Z] = FunR3.XYZ;
PW = exp(1i*k*(u(1)*X + u(2)*Y + u(3)*Z));
RHS = integral(Gamma,Vh,PW);

if is2d(m)
    if k==0
        Gxy = @(X,Y)(femGreenKernel(X,Y,'[log(r)]',0));
        dGxy = cell(3,1);
        dGxy{1} = @(X,Y)(femGreenKernel(X,Y,'grady[log(r)]1',0));
        dGxy{2} = @(X,Y)(femGreenKernel(X,Y,'grady[log(r)]2',0));
        dGxy{3} = @(X,Y)(femGreenKernel(X,Y,'grady[log(r)]3',0));
        
        C = -1/(2*pi);
    else
        Gxy = @(X,Y)(femGreenKernel(X,Y,'[H0(kr)]',k));
        dGxy = cell(3,1);
        dGxy{1} = @(X,Y)(femGreenKernel(X,Y,'grady[H0(kr)]1',k));
        dGxy{2} = @(X,Y)(femGreenKernel(X,Y,'grady[H0(kr)]2',k));
        dGxy{3} = @(X,Y)(femGreenKernel(X,Y,'grady[H0(kr)]3',k));
        C = 1i/4;
    end
    regKer = '[log(r)]';
    dregKer = 'grady[log(r)]';
    Creg = -1/(2*pi);
else
    if k==0
        Gxy = @(X,Y)(femGreenKernel(X,Y,'[1/r]',0));
        dGxy = cell(3,1);
        dGxy{1} = @(X,Y)(femGreenKernel(X,Y,'grady[1/r]1',0));
        dGxy{2} = @(X,Y)(femGreenKernel(X,Y,'grady[1/r]2',0));
        dGxy{3} = @(X,Y)(femGreenKernel(X,Y,'grady[1/r]3',0));
    else
        Gxy = @(X,Y)(femGreenKernel(X,Y,'[exp(ikr)/r]',k));
        dGxy = cell(3,1);
        dGxy{1} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k));
        dGxy{2} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k));
        dGxy{3} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k));
        
    end
    regKer = '[1/r]';
    dregKer = 'grady[1/r]';
    C = 1/(4*pi);
    Creg = C;
end




Sk = C*integral(Gamma,Gamma,Vh,Gxy,Vh) ...
    + Creg*regularize(Gamma,Gamma,Vh,regKer,Vh);

Dk = C*integral(Gamma,Gamma,Vh,dGxy,ntimes(Vh))...
    + Creg*regularize(Gamma,Gamma,Vh,dregKer,ntimes(Vh));

Mass = integral(Gamma,Vh,Vh);

eta = -1i*k/2;

A = @(x)(eta*Sk*x + (Mass*x/2 + Dk*x));

t = tic;
[lambda,flag,~,iter,resvec] = gmres(A,RHS,[],1e-8,size(Sk,1));
tSolve = toc(t);
resvec = resvec/norm(RHS);

niter = iter(end);
if flag ~=0
    error('Gmres did not converge');
end


end

