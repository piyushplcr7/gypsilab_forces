clear; clc;
t=tic;
mesh = mshSphere(200,1);
S1_Gamma = fem(mesh,'P1');
kernelx = @(x,y,z) 1./(vecnorm(z,2,2));
t2mat = panel_oriented_assembly(mesh,kernelx,curl(S1_Gamma),curl(S1_Gamma));
elt=toc(t);