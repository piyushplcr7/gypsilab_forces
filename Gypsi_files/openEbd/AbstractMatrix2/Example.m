%% Example of use
% Create an Abstract matrix and solve the associated linear system
clear all
close all
clc;

targetNnz = 5*10^2;
m = 10^3; density = targetNnz/(m*m);
concretePart = sprand(m,m,density);

lambda = rand(m,1);
abstractPart = @(x)(exp(lambda)*sum(exp(-lambda).*x));


M = AbstractMatrix(concretePart,abstractPart);

M = M + 25*AbstractMatrix.eye(m) + 0.01*AbstractMatrix.ones(m) +...
    AbstractMatrix.spdiag(randn(m,1));
% eye and ones are stored in a pure abstract way by default. 

spy(M);
title('Sparsity pattern of the concrete part of M');
drawnow;

x0 = randn(m,1);
x0 = x0/norm(x0);
b = M*x0;
[x1,flag,relres,iter,resvec] = gmres(M,b); % Uses the concrete part as a preconditioner by default. 

figure
semilogy(1:length(resvec),resvec,'-o');
xlabel('Iteration number')
ylabel('Relative residual')
title('Gmres convergence hisotry')

relerr1 = norm(x1-x0)/norm(x0); 
fprintf('Solution found with relative error %s \n\n',num2str(relerr1));

Mf = full(M);
x2 = Mf\b;

relerr2 = norm(x2-x0)/norm(x0); 
fprintf('Relative error using full matrix %s \n\n',num2str(relerr2));

imagesc(Mf);
title('Representation of the full Matrix')