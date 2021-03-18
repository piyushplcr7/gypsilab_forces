%% Represent spherical harmonics


clean; 

m = mshSphere(2000,1);
Vh = P1(m);
Gamma = dom(m,7);

tic;

M = integral(Gamma,Vh,Vh);
K = integral(Gamma,grad(Vh),grad(Vh));

toc;

[P,D] = eig(full(M\K));
[d,I] = sort(diag(D),'ascend');
d(1:10)
N = 6;
for i = 1:N
    for j = 1:N
        k = (i-1)*N + j;
        subplot(N,N,k)
        
        vals = feval(Vh,P(:,I(k)),m);
        vals = vals./max(abs(vals));
        s = vals;
        s(abs(vals) < 1e-5) = 0; 
        rho = abs(vals);
        m2 = m; 
        m2.vtx = rho.*m.vtx./sqrt((scal3D(m.vtx,m.vtx)));
        plot(m2,s);
        axis equal
        view(49,25);
        caxis([-1,1])
        title(sprintf('lambda = %s',num2str(d(k))));
    end
end
    