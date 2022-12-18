% Double layer test

clear;clc;

mesh = mshSphere(10,1);
Gamma = dom(mesh,3);
NED = fem(mesh,'NED');
DIV = fem(mesh,'RWG');
P1 = fem(mesh,'P1');
DIV0 = nxgrad(P1);

C = double_layer(Gamma,DIV,NED);

kernel = cell(3,1);
kernel{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
kernel{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
kernel{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);

MK = 1/4/pi*integral(Gamma,Gamma,DIV,kernel,DIV);
MK = MK + 1/4/pi*regularize(Gamma,Gamma,DIV,'grady[1/r]',DIV);

norm(MK+C)



