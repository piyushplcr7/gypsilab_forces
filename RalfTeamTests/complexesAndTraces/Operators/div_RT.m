function [divU,Wh] = div_RT(Vh)

% Returns the P0 element equal to div u.
m = Vh.msh;
Wh = fem(m,'P0');

gss    = 1;
Omega = dom(m,gss);
M = integral(Omega,Wh,Wh);
Fv = integral(Omega,Wh,div(Vh));

divU = spdiags(1./diag(M),0,size(M,1),size(M,2))*Fv;


end
