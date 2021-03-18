% Poisson equation with Raviart-Thomas elements.

clean;

m = mshSquare(500,[1,1]);
Gamma = dom(m,3);
RT = fem(m,'RWG');
P0 = fem(m,'P0');
P1 = fem(m,'P1');
P10 = dirichlet(fem(m,'P1'),m.bnd);


z = @(X)(0*X(:,1));
f = @(X)(sin(5*X(:,1).*X(:,2)));

M11 = integral(Gamma,RT,RT);
M12 = integral(Gamma,P0,div(RT)).';
M21 = -integral(Gamma,div(RT),P0).';
M22 = integral(Gamma,P0,P0);


N1 = size(M21,1);
N2 = size(M12,2);

M = [M11 M12; M21 M22];

U1 = zeros(size(M11,1),1);
U2 = integral(Gamma,P0,f);
U = [U1;U2];

n1 = length(U1);
n2 = length(U2);

sol = M\U;
u = sol((n1+1):end);


m2 = m;
m2.vtx(:,3) = feval(P0,u,m);
plot(m2,m2.vtx(:,3));


%% Comparison with classical formulation

M = integral(Gamma,grad(P1),grad(P1)) + integral(Gamma,P1,P1);
U = integral(Gamma,P1,f);

sol = M\U;

m2 = m;
m2.vtx(:,3) = feval(P1,sol,m);
figure;
plot(m2,m2.vtx(:,3));