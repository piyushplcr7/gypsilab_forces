function [lambda,niter,resvec,tSolve] = solveSingleLayer(m,k,u,varargin)

p = inputParser;
p.addOptional('niter',size(m.vtx,1));
p.addOptional('precond','none');
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

switch p.Results.precond
    case 'none'
        Prec = @(x)(x);
    case 'Calderon'
        Nk = C*integral(Gamma,Gamma,nxgrad(Vh),Gxy,nxgrad(Vh)) ...
            + Creg*regularize(Gamma,Gamma,nxgrad(Vh),regKer,nxgrad(Vh));
        Nk = Nk - k^2*(C*integral(Gamma,Gamma,ntimes(Vh),Gxy,ntimes(Vh))...
            +Creg*regularize(Gamma,Gamma,ntimes(Vh),regKer,ntimes(Vh)));
        Mass = integral(Gamma,Vh,Vh);
        [L,U] = lu(Mass);
        invM = @(x)(U\(L\x));
        Prec = @(x)(invM(Nk*invM(x)));
    case 'S0'
        Gxy0 = @(X,Y)(femGreenKernel(X,Y,regKer,0));
        S0 = Creg*integral(Gamma,Gamma,Vh,Gxy0,Vh) ...
            + Creg*regularize(Gamma,Gamma,Vh,regKer,Vh);
        Prec = @(x)(S0\x);
        
    case 'CFIE'
        
        
        Sk = C*integral(Gamma,Gamma,Vh,Gxy,Vh) ...
            + Creg*regularize(Gamma,Gamma,Vh,regKer,Vh);
        
        Dk = C*integral(Gamma,Gamma,Vh,dGxy,ntimes(Vh))...
            + Creg*regularize(Gamma,Gamma,Vh,dregKer,ntimes(Vh));
        
        Mass = integral(Gamma,Vh,Vh);
        eta = 1i*k/2;
        
        CFIE = @(x)(eta*Sk*x - (Mass*x/2 + Dk*x));
        t = tic;
        [lambda,flag,relres,niter,resvec] = gmres(CFIE,RHS,[],1e-8,size(Sk,1));
        niter = niter(end);
        tSolve = toc(t);
        radiat = mshSquare(500,[4,4]);
        SL = C*integral(radiat.vtx,Gamma,Gxy,Vh)...
            + Creg*regularize(radiat.vtx,Gamma,regKer,Vh);
        
        DL = C*integral(radiat.vtx,Gamma,dGxy,ntimes(Vh))...
            + Creg*regularize(radiat.vtx,Gamma,dregKer,ntimes(Vh));
        
        rad = SL*(eta*lambda) - DL*(lambda);
        figure;
        plotOn(radiat,real(rad));
        hold on
        plot(m);
        return
    case 'Darbas'
        Delta = integral(Gamma,grad(Vh),grad(Vh));
        Mass = integral(Gamma,Vh,Vh);
        Dk = C*integral(Gamma,Gamma,Vh,dGxy,ntimes(Vh))...
            + Creg*regularize(Gamma,Gamma,Vh,dregKer,ntimes(Vh));
        I = eye(size(Dk,1));
        X = m.vtx(1,:);
        r = max(norm3D(X - m.vtx));
        keps = k + 1i*0.4*r*k^(1/3);
        [L,U] = lu(Mass);
        invM = @(x)(U\(L\x));
        PrecDarbas = @(x)(padePrecondDarbas(x,15,pi/3,keps,Mass,Delta));
        GCSIE = @(x)(invM(Sk*invM(PrecDarbas(x))) - I/2*x - invM(Dk*x));
        
        t = tic;
        [lambda,flag,relres,niter,resvec] = gmres(GCSIE,invM(RHS),[],1e-8,size(Sk,1));
        tSolve = toc(t);
        niter = niter(end);
        
        radiat = mshSquare(500,[4,4]);
        SL = C*integral(radiat.vtx,Gamma,Gxy,Vh)...
            + Creg*regularize(radiat.vtx,Gamma,regKer,Vh);
        
        DL = C*integral(radiat.vtx,Gamma,dGxy,ntimes(Vh))...
            + Creg*regularize(radiat.vtx,Gamma,dregKer,ntimes(Vh));
        
        rad = SL*(invM(PrecDarbas(lambda))) - DL*(lambda);
        figure;
        plotOn(radiat,real(rad));
        hold on
        plot(m);
        
        
        return
end


t = tic;
[lambda,flag,relres,iter,resvec] = gmres(@(x)(Prec(Sk*x)),Prec(RHS),[],1e-8,size(Sk,1));
tSolve = toc(t);

radiat = mshSquare(500,[4,4]);
SL = C*integral(radiat.vtx,Gamma,Gxy,Vh)...
    + Creg*regularize(radiat.vtx,Gamma,regKer,Vh);

radiat2 = rotate(mshSquare(500,[4,4]),[1 0 0],pi/2);
radiat2 = translate(radiat2,[0 2 0]);

SL2 = C*integral(radiat2.vtx,Gamma,Gxy,Vh)...
    + Creg*regularize(radiat2.vtx,Gamma,regKer,Vh);

rad = SL*(lambda);
rad2 = SL2*(lambda);
figure;
plotOn(radiat,real(rad));
hold on
plotOn(radiat2,real(rad2));
plot(m);
plotOn(m,real(lambda));



resvec = resvec./norm(Prec(RHS));

niter = iter(end);
if flag ~=0
    error('Gmres did not converge');
end


end

