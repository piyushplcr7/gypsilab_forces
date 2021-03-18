function M = infExt(fe)
m = fe.msh;
dm = mshBoundary(m);
fe0 = dirichlet(fe,dm);
dfe = fem(dm,fe.typ);
gss    = 3;
Gamma = dom(dm,gss);
bool   = (size(m.elt,2) == 4);
gss2    = 4*bool + 3*~bool;
Omega = dom(m,gss2);
Mtmp = trace(fe).';
a = integral(Omega,grad(fe0),grad(fe0));
gradfe1 = grad(fe0);
gradfe2 = grad(fe);
Uqm1 = gradfe1.uqm(Omega);
Uqm2 = gradfe2.uqm(Omega);
[X,W] = Omega.qud;
W = spdiags(W,0,size(W,1),size(W,1));
b = sparse([],[],[],size(Uqm1{1},2),size(Mtmp,2));
for i = 1:3
    b= b - Uqm1{i}.'*W*Uqm2{i}*Mtmp;
end
ExtHomogeneous = a\b;
M = Mtmp + (integral(Omega,fe,fe)\integral(Omega,fe,fe0))*ExtHomogeneous;
end