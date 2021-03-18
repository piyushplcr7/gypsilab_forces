function [divnxU,Wh] = divnx_NED(Vh)

assert(size(Vh.msh.elt,2)==3);
Wh = fem(Vh.msh,'P0');
Gamma = dom(Vh.msh,1);
M = integral(Gamma,Wh,Wh);
Fv = integral(Gamma,Wh,curl(Vh));
divnxU = spdiags(1./diag(M),0,size(M,1),size(M,2))*Fv;

end