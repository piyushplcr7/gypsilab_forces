% Benchmark Cube
close all;
clear all;

Ns = [25 50 100 200 400 600 800];




[X,Y,Z] = FunR3.XYZ;


for i = 1:length(Ns)
    
    N = Ns(i);
    
    
    m = bnd(mshCube(N,[1 1 1]));
    m = rotate(m,[1 1 1],pi/4);
    wall = rotate(mshSquare(800,[4,4]),[1 0 0],pi/2);
    wall = translate(wall,[0 2 0]);
    
    stp = m.stp;
    hh = stp(1);
    k = 1.2/hh
    Gk = @(X,Y)(femGreenKernel(X,Y,'[exp(ikr)/r]',k));
    G0 = @(X,Y)(femGreenKernel(X,Y,'[1/r]',0));
    dGk = cell(3,1);
    dGk{1} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k));
    dGk{2} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k));
    dGk{3} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k));
    
    PW = exp(1i*k*(Y + Z)/sqrt(2));
    
    Vh = fem(m,'P1');
    Gamma = dom(m,3);
    
    Sk = 1/(4*pi)*integral(Gamma,Gamma,Vh,Gk,Vh) ...
        + 1/(4*pi)*regularize(Gamma,Gamma,Vh,'[1/r]',Vh);
    
    Dk = 1/(4*pi)*integral(Gamma,Gamma,Vh,dGk,ntimes(Vh))...
        + 1/(4*pi)*regularize(Gamma,Gamma,Vh,'grady[1/r]',ntimes(Vh));
    
    Nk = 1/(4*pi)*integral(Gamma,Gamma,nxgrad(Vh),Gk,nxgrad(Vh)) ...
        + 1/(4*pi)*regularize(Gamma,Gamma,nxgrad(Vh),'[1/r]',nxgrad(Vh))...
        -k^2/(4*pi)*integral(Gamma,Gamma,ntimes(Vh),Gk,ntimes(Vh))...
        -k^2/(4*pi)*regularize(Gamma,Gamma,ntimes(Vh),'[1/r]',ntimes(Vh));
    
    N0 = 1/(4*pi)*integral(Gamma,Gamma,nxgrad(Vh),G0,nxgrad(Vh)) ...
        + 1/(4*pi)*regularize(Gamma,Gamma,nxgrad(Vh),'[1/r]',nxgrad(Vh));
    
    
    M = integral(Gamma,Vh,Vh);
    Delta = integral(Gamma,grad(Vh),grad(Vh));
    
    keps = k + 1i*0.5*k^(1/3);
    DtNk = -padePrecondDarbas([],15,pi/2,keps,M,Delta);
    I = eye(size(M,1));
    CFIE = -(M\Sk)*(1i*k) + (I/2 + M\Dk);
    GCSIE = (M\Sk)*(M\DtNk) + (I/2 + M\Dk);
    
    NoPrec = eye(size(Sk,1));
    Prec0 = M^(-1)*(N0+M)*M^(-1);
    Preck = M^(-1)*(Nk)*M^(-1);
    
    PW = exp(1i*k*(Y+Z)/sqrt(2));
    
    rhs = integral(Gamma,Vh,PW);
    
    maxit = max(500,size(Sk,1));
    
    [u,~,~,~,RESVEC1] = gmres(NoPrec*Sk,NoPrec*rhs,[],1e-8,maxit);
    [~,~,~,~,RESVEC2] = gmres(Prec0*Sk,Prec0*rhs,[],1e-8,maxit);
    [~,~,~,~,RESVEC3] = gmres(Preck*Sk,Preck*rhs,[],1e-8,maxit);
    [~,~,~,~,RESVEC4] = gmres(CFIE,M\rhs,[],1e-8,maxit);
    [~,~,~,~,RESVEC5] = gmres(GCSIE,M\rhs,[],1e-8,maxit);
    
    iterNoPrec(i) = size(RESVEC1,1)
    iterN0Prec(i) = size(RESVEC2,1)
    iterNkPrec(i) = size(RESVEC3,1)
    iterCFIEPrec(i) = size(RESVEC4,1)
    iterGCSIEPrec(i) = size(RESVEC5,1)
    
    stp = m.stp;
    h(i) = stp(3);
    
end


figure;
loglog(h,iterNoPrec,'-x')
hold on
loglog(h,iterN0Prec,'-x')
loglog(h,iterNkPrec,'-x')
loglog(h,iterCFIEPrec,'-x')
loglog(h,iterGCSIEPrec,'-x')
grid on
xlabel('$h$','Interpreter','latex')
ylabel('$n_{it}(1$e-$8)$','Interpreter','latex')
legend({'No Prec.','N0','Nk','CFIE','GCSIE'})
ylim([min(iterGCSIEPrec)-3,max(iterNoPrec)+3])


