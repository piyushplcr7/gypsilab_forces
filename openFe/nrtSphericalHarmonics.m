
clean; 

m = mshSphere(200,1);
Vh = P2(m);
Gamma = dom(m,7);

tic;

M = integral(Gamma,Vh,Vh);
K = integral(Gamma,grad(Vh),grad(Vh));

toc;

[P,D] = eig(full(M\K));
[d,I] = sort(diag(D),'ascend');
d(1:10)
N = 5;
for i = 1:N
    for j = 1:N
        k = (i-1)*N + j;
        subplot(N,N,k)
        vals = feval(Vh,P(:,I(k)),m);
        plot(m,vals);
    end
end
    

Vh = fem(m,'P2');
Gamma = dom(m,7);

tic;

M = integral(Gamma,Vh,Vh);
K = integral(Gamma,grad(Vh),grad(Vh));

toc;

[P,D] = eig(full(M\K));
[d,I] = sort(diag(D),'ascend');
d(1:10)
figure;
N = 5;
for i = 1:N
    for j = 1:N
        k = (i-1)*N + j;
        subplot(N,N,k)
        vals = feval(Vh,P(:,I(k)),m);
        plot(m,vals);
    end
end
    