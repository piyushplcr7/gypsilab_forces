close all;
clear all;
clc;

%% Frequency fixed

Ns = [100 200 400 800 1600 3200];

k = 5;

for n = 1:length(Ns)
    tol = 1e-3;
    N = Ns(n);
    m = mshTorus(N,3,1);
    
    Gxy = @(X,Y)(femGreenKernel(X,Y,'[exp(ikr)/r]',k));
    Gamma = dom(m,3);
    Vh = fem(m,'P1');
    
    Sk = 1/(4*pi)*(integral(Gamma,Gamma,Vh,Gxy,Vh,tol) ...
        + regularize(Gamma,Gamma,Vh,'[1/r]',Vh));
    
    
    [X,Y,Z] = FunR3.XYZ;
    u = [1 1 1];
    u = u./norm(u,2);
    incWave = exp(1i*k*(u(1)*X + u(2)*Y + u(3)*Z));
    
    L = integral(Gamma,Vh,incWave);
    
    
    [lambda,FLAG,RELRES,ITER,RESVEC] = gmres(@(x)(Sk*x),L,[],1e-8,size(Sk,1));
    
    niter(n) = ITER(end);
    
    
end


figure;
plotOn(m,real(lambda));

plot(Ns,niter,'r*')
hold on
plot(Ns,niter)


%% N fixed

close all;
clear all;
clc;


N = 1600;



m = mshTorus(N,3,1);
Gamma = dom(m,3);
Vh = fem(m,'P1');
reg = regularize(Gamma,Gamma,Vh,'[1/r]',Vh);


ks = [1 2 3 4 5 6 7 8];

for n = 1:length(ks)
    tol = 1e-3;
    k = ks(n);
    
    Gxy = @(X,Y)(femGreenKernel(X,Y,'[exp(ikr)/r]',k));
    
    
    Sk = 1/(4*pi)*(integral(Gamma,Gamma,Vh,Gxy,Vh,tol) ...
        + reg);
    
    
    [X,Y,Z] = FunR3.XYZ;
    u = [1 1 1];
    u = u./norm(u,2);
    incWave = exp(1i*k*(u(1)*X + u(2)*Y + u(3)*Z));
    
    L = integral(Gamma,Vh,incWave);
    
    
    [lambda,FLAG,RELRES,ITER,RESVEC] = gmres(@(x)(Sk*x),L,[],1e-8,size(Sk,1));
    
    niter(n) = ITER(end);
    
    
end

figure;
plot(ks,niter,'r*')
hold on
plot(ks,niter)
xlabel('Ndof')
ylabel('niter GMRES')


%% Resolve frequency

clear all;
close all;
clc;

ks = [1 2 3 4 6 8 16]/2;
Ns = fix((2*pi)^2*ks.^2);




for n = 1:length(ks)
    tol = 1e-3;
    k = ks(n);
    N = Ns(n);
    m = mshTorus(N,3,1);
    Gamma = dom(m,3);
    Vh = fem(m,'P1');
    Gxy = @(X,Y)(femGreenKernel(X,Y,'[exp(ikr)/r]',k));
    
    
    Sk = 1/(4*pi)*(integral(Gamma,Gamma,Vh,Gxy,Vh) ...
        + regularize(Gamma,Gamma,Vh,'[1/r]',Vh));
    
    
    [X,Y,Z] = FunR3.XYZ;
    u = [1 1 1];
    u = u./norm(u,2);
    incWave = exp(1i*k*(u(1)*X + u(2)*Y + u(3)*Z));
    
    L = integral(Gamma,Vh,incWave);
    
    tmp = tic;
    [lambda,FLAG,RELRES,ITER,RESVEC] = gmres(@(x)(Sk*x),L,[],1e-8,size(Sk,1));
    t(n) = toc(tmp);
    
    niter(n) = ITER(end);
    
    
end

figure;
plotOn(m,real(lambda));
axis equal
axis off;

figure
loglog(Ns,niter);
hold on
loglog(Ns,niter,'r*');
xlabel('Ndof')
ylabel('niter GMRES')


figure;
plot(Ns,t);
hold on
plot(Ns,t,'r*');
xlabel('Ndof')
ylabel('Resolution time(s)')


%% Calderon preconditioning


clear all;
close all;
clc;

ks = [1 2 3 4 6 8]/2;
Ns = fix((2*pi)^2*ks.^2);




for n = 1:length(ks)
    tol = 1e-3;
    k = ks(n);
    N = Ns(n);
    m = mshSphere(N,1);
    Gamma = dom(m,3);
    Vh = fem(m,'P1');
    Gxy = @(X,Y)(femGreenKernel(X,Y,'[exp(ikr)/r]',k));
    
    
    Sk = 1/(4*pi)*(integral(Gamma,Gamma,Vh,Gxy,Vh) ...
        + regularize(Gamma,Gamma,Vh,'[1/r]',Vh));
    
    Nk = 1/(4*pi)*(integral(Gamma,Gamma,nxgrad(Vh),Gxy,nxgrad(Vh))...
        - k^2*integral(Gamma,Gamma,ntimes(Vh),Gxy,ntimes(Vh))...
        + regularize(Gamma,Gamma,nxgrad(Vh),'[1/r]',nxgrad(Vh))...
         - k^2*regularize(Gamma,Gamma,ntimes(Vh),'[1/r]',ntimes(Vh)));
    
    Mass = integral(Gamma,Vh,Vh);
    [L,U] = lu(Mass);
    invM = @(x)(U\(L\x));
    
    [X,Y,Z] = FunR3.XYZ;
    u = [1 1 1];
    u = u./norm(u,2);
    incWave = exp(1i*k*(u(1)*X + u(2)*Y + u(3)*Z));
    
    L = integral(Gamma,Vh,incWave);
    
    tmp = tic;
    [lambda,FLAG,RELRES,ITER,RESVEC] = gmres(@(x)(invM(Nk*invM(Sk*x))),invM(Nk*invM(L)),[],1e-8,size(Sk,1));
    t(n) = toc(tmp);
    
    niter(n) = ITER(end);
    
    
end

figure;
plotOn(m,real(lambda));
axis equal
axis off;

figure
loglog(Ns,niter);
hold on
loglog(Ns,niter,'r*');
xlabel('Ndof')
ylabel('niter GMRES')


figure;
plot(Ns,t);
hold on
plot(Ns,t,'r*');
xlabel('Ndof')
ylabel('Resolution time(s)')



