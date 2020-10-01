%% Eigenvalues of laplace operator on the disk.

%% New fe

clean;

m = mshDisk(800,1);
dir = m.bnd;
Vh = dirichlet(P1(m),dir);
Gamma = dom(m,3);


I = integral(Gamma,Vh,Vh);
K = integral(Gamma,grad(Vh),grad(Vh));


[P,D] = eig(full(I\K));
[d,I] = sort(diag(D),'ascend');
d(1:10)

vals = feval(Vh,P(:,I(7)),m);

figure;
m2 = m;
m2.vtx(:,3) = vals;
plot(m2,vals);
view(45,45);

%% Old fem

dir = m.bnd;
Wh = dirichlet(fem(m,'P1'),dir);
Gamma = dom(m,3);

I = integral(Gamma,Wh,Wh);
K = integral(Gamma,grad(Wh),grad(Wh));


[P,D] = eig(full(I\K));
[d,I] = sort(diag(D),'ascend');
d(1:10)
vals = feval(Wh,P(:,I(7)),m);

figure;
m2 = m;
m2.vtx(:,3) = vals;
plot(m2,vals);
view(45,45);

