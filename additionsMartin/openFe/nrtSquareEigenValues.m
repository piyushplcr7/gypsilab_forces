% 
clean;

m = mshSquare(500,[pi,pi]);


Vh = P2(m);
Gamma = dom(m,7);

M1 = integral(Gamma,Vh,Vh);
K1 = integral(Gamma,grad(Vh),grad(Vh));


[P,D] = eig(full(M1\K1));
[d,I] = sort(diag(D),'ascend');
d(1:10)
vals = real(feval(Vh,P(:,I(12)),m));
m2 = m;
m2.vtx(:,3) = vals;
plot(m2,vals);
view(45,45);
err1 = norm(d(1:10) - [0;1;1;2;4;4;5;5;8;9],2);
assert(err1<0.01);




%% Comparison 


Vh = fem(m,'P2');
Gamma = dom(m,7);

M2 = integral(Gamma,Vh,Vh);
K2 = integral(Gamma,grad(Vh),grad(Vh));


[P,D] = eig(full(M2\K2));
[d,I] = sort(diag(D),'ascend');
d(1:10)

vals = real(feval(Vh,P(:,I(12)),m));
m.vtx(:,3) = vals;
figure;
plot(m,vals);
view(45,45);
err2 = norm(d(1:10) - [0;1;1;2;4;4;5;5;8;9],2);
assert(err2<0.01);
