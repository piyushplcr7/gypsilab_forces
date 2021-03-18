function M = trace(fe)
% if fe.typ is Pk
m = fe.msh;
dm = mshBoundary(m);
dfe = fem(dm,fe.typ);
gss    = 3;
Gamma = dom(dm,gss);
M1 = integral(Gamma,dfe,dfe);
M2 = integral(Gamma,dfe,fe);
M = M1\M2;
end

