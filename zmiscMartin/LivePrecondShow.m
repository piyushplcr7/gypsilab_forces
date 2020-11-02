%% Live precond show



%% Sphere, fixed k, increasing N. 

clear all
close all
clc;
k = 4;


Ns = [10 20 40 80 160];
u = [1 1 1]/sqrt(3);
figure;
resvec = cell(length(Ns),1);
for n = 1:length(Ns)
    m = mshSphere(Ns(n),1);
    [lambda,niter(n),resvec{n}] = solveSingleLayer(m,k,u);
    resvec{n}(resvec{n}==0) = 1e-15;
    semilogy(1:(niter(n)+1),resvec{n},'DisplayName',num2str(Ns(n)));
    hold on
end

xlabel('iteration count')
ylabel('residual error')
legend show
title('Gmres convergence comparison')

figure
plot(Ns,niter);
xlabel('Ndof'); ylabel('Gmres it');
title('Number of GMRES iteration')

figure;
plotOn(m,real(lambda));
title('Sample solution')

%% Sphere increasing k and N

clear all;
close all;
clc

ks = [1 2 4 8 16]/4;

Ns = [10 20 40 80 160];
u = [1 1 1]/sqrt(3);
figure;
resvec = cell(length(Ns),1);
for n = 1:length(Ns)
    k = ks(n);
    m = mshSphere(Ns(n),1);
    [lambda,niter(n),resvec{n}] = solveSingleLayer(m,k,u);
    resvec{n}(resvec{n}==0) = 1e-15;
    semilogy(1:(niter(n)+1),resvec{n},'DisplayName',num2str(Ns(n)));
    hold on
end
xlabel('iteration count')
ylabel('residual error')
legend show
title('Gmres convergence comparison')

figure
plot(Ns,niter);
xlabel('Ndof'); ylabel('Gmres it');
title('Number of GMRES iteration')

figure;
plotOn(m,real(lambda));
title('Sample solution')

%% Sphere, increasing k and N, Calderon prec

clear all;
close all;
clc;


ks = [1 2 4 8 16]/4;

Ns = [10 20 40 80 160];
u = [1 1 1]/sqrt(3);
figure;
resvec = cell(length(Ns),1);
for n = 1:length(Ns)
    k = ks(n);
    m = mshSphere(Ns(n),1);
    [lambda,niter(n),resvec{n}] = solveSingleLayer(m,k,u,...
        'precond','Calderon');
    resvec{n}(resvec{n}==0) = 1e-15;
    semilogy(1:(niter(n)+1),resvec{n},'DisplayName',num2str(Ns(n)));
    hold on
end

xlabel('iteration count')
ylabel('residual error')
legend show
title('Gmres convergence comparison')

figure
plot(Ns,niter);
xlabel('Ndof'); ylabel('Gmres it');
title('Number of GMRES iteration')

figure;
plotOn(m,real(lambda));
title('Sample solution')

%% Torus, increasing k and N, Calderon prec

clear all;
close all;
clc;

ks = [2 4 8 16 32]/2;
Ns = max(2*ks.^2,20);
u = [1 1 1]/sqrt(3);
resvec = cell(length(Ns),1);
for n = 1:length(Ns)
    k = ks(n);
    m = mshSphere(Ns(n),1);
%     m = m.rotate([1 2 -3],pi/6);
    [lambda,niter(n),resvec{n},tSolve] = solveSingleLayer(m,k,u,...
        'precond','Calderon');
    [lambda,niter2(n),resvec{n},tSolve] = solveSingleLayer(m,k,u,...
        'precond','S0');
    [lambda,niter3(n),resvec2{n},tSolve2] = solveSingleLayer(m,k,u,...
        'precond','none');
end

h = figure;
loglog(ks,niter,'DisplayName','Calderon');
hold on
loglog(ks,niter2,'DisplayName','S_0^{-1}');
loglog(ks,niter3,'DisplayName','No precond.');
loglog(ks,5*ks,'k--','DisplayName','O(k)');

xlabel('k'); ylabel('Gmres it');
title('Number of GMRES iteration')
legend show

figure;
plotOn(m,real(lambda));
title('Sample solution')
axis equal;


%% Torus, add pseud-diff preconds
% clear all
% close all
% clc;

ks = [2 4 8 16 32]/2;
Ns = max(2*ks.^2,20);

u = [1 1 1]/sqrt(3);
niter1 = zeros(length(Ns),1);
niter2 = zeros(length(Ns),1);
for n = 1:length(Ns)
    k = ks(n);
    m = mshSphere(Ns(n),1);
    [lambda,niter1(n),resvec,tSolve] = solveSingleLayer(m,k,u,...
        'precond','CFIE');
    [lambda,niter2(n),resvec,tSolve] = solveSingleLayer(m,k,u,...
        'precond','Darbas');
end

figure(h);
hold on
loglog(ks,niter1,'DisplayName','CFIE');
loglog(ks,niter2,'DisplayName','GCSIE');



%% At least, Calderon is an optimal preconditioning method

clear all;
close all;
clc;

ks = [2 2 2 2 2];
Ns = [40 80 160 320 640];
u = [1 1 1]/sqrt(3);
figure;
resvec = cell(length(Ns),1);
for n = 1:length(Ns)
    k = ks(n);
    m = mshTorus(Ns(n),1,0.3);
    [lambda,niter(n),resvec{n},tSolve] = solveSingleLayer(m,k,u,...
        'precond','Calderon');
    resvec{n}(resvec{n}==0) = 1e-15;
    semilogy(1:(niter(n)+1),resvec{n},'DisplayName',num2str(Ns(n)));
    hold on
end

xlabel('iteration count')
ylabel('residual error')
legend show
title('Gmres convergence comparison')

figure
plot(Ns,niter);
xlabel('Ndof'); ylabel('Gmres it');
title('Number of GMRES iteration')

figure;
plotOn(m,real(lambda));
title('Sample solution')
axis equal;


 
%% A convex domain then 

clear all;
close all;
clc;
ks = [2 4 8 16 32 64]/3;

Ns = [40 80 160 320 640];
u = [0.5 2 1]/sqrt(5.25);
figure;
resvec = cell(length(Ns),1);
for n = 1:length(Ns)
    k = ks(n);
    m = mshCube(Ns(n),[1 1 1]);
    m = m.bnd;
    [lambda,niter(n),resvec{n},tSolve] = solveSingleLayer(m,k,u,...
        'precond','Calderon');
    resvec{n}(resvec{n}==0) = 1e-15;
    semilogy(1:(niter(n)+1),resvec{n},'DisplayName',num2str(Ns(n)));
    hold on
end

xlabel('iteration count')
ylabel('residual error')
legend show
title('Gmres convergence comparison')

figure
plot(Ns,niter);
xlabel('Ndof'); ylabel('Gmres it');
title('Number of GMRES iteration')

figure;
plotOn(m,real(lambda));
title('Sample solution')
axis equal;
 
%% Cube, last example, no preconditioner

clc
[~,niterLast,resvec{n},tSolveNoPrec] = solveSingleLayer(m,k,u,...
        'precond','none');
    
 fprintf('%s it (no prec) vs %s it (prec) \n',num2str(niterLast),num2str(niter(end)));
fprintf('%s s (no prec) vs %s s (prec) \n',num2str(tSolve),num2str(tSolveNoPrec));
 
 %% CFIE for the torus
 
clear all;
close all;
clc;

ks = [1 2 4 8 16]/4;

Ns = [40 80 160 320 640];
u = [1 1 1]/sqrt(3);
figure;
resvec = cell(length(Ns),1);
for n = 1:length(Ns)
    k = ks(n);
    m = mshTorus(Ns(n),1,0.3);
    [~,niter(n),resvec{n}] = solveCFIE(m,k,u);
    resvec{n}(resvec{n}==0) = 1e-15;
    semilogy(1:(niter(n)+1),resvec{n},'DisplayName',num2str(Ns(n)));
    hold on
end
xlabel('iteration count')
ylabel('residual error')
legend show
title('Gmres convergence comparison')

figure
plot(Ns,niter);
xlabel('Ndof'); ylabel('Gmres it');
title('Number of GMRES iteration')


%% Is CFIE an optimal preconditioning strategy ?

clear all;
close all;
clc;

ks = [4 4 4 4 4];

Ns = [40 80 160 320 640];
u = [1 1 1]/sqrt(3);
figure;
resvec = cell(length(Ns),1);
for n = 1:length(Ns)
    k = ks(n);
    m = mshTorus(Ns(n),1,0.3);
    [~,niter(n),resvec{n}] = solveCFIE(m,k,u);
    resvec{n}(resvec{n}==0) = 1e-15;
    semilogy(1:(niter(n)+1),resvec{n},'DisplayName',num2str(Ns(n)));
    hold on
end
xlabel('iteration count')
ylabel('residual error')
legend show
title('Gmres convergence comparison')

figure
plot(Ns,niter);
xlabel('Ndof'); ylabel('Gmres it');
title('Number of GMRES iteration')




%% Does CFIE work on a cube ?

clear all;
close all;
clc;

ks = [1 2 4 8 16]/4;

Ns = ks.^2*50;
u = [1 1 1]/sqrt(3);
figure;
resvec = cell(length(Ns),1);
for n = 1:length(Ns)
    k = ks(n);
    m = bnd(mshCube(Ns(n),[1 1 1]));
    [~,niter(n),resvec{n}] = solveCFIE(m,k,u);
    resvec{n}(resvec{n}==0) = 1e-15;
    semilogy(1:(niter(n)+1),resvec{n},'DisplayName',num2str(Ns(n)));
    hold on
end
xlabel('iteration count')
ylabel('residual error')
legend show
title('Gmres convergence comparison')

figure
plot(Ns,niter);
xlabel('Ndof'); ylabel('Gmres it');
title('Number of GMRES iteration')



 
 