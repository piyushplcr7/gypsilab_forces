clear all
close all


m = mshSphere(500,1);mplot = mshSquare(500,[1,1]);X = mplot.vtx;

k = 15;
G = @(X,Y)(femGreenKernel(X,Y,'[exp(ikr)/r]',k));
dG{1} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]1',k));
dG{2} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]2',k));
dG{3} = @(X,Y)(femGreenKernel(X,Y,'grady[exp(ikr)/r]3',k));
Vh = fem(m,'P1');
Gamma = dom(m,3);
SL = 1i/4*(integral(X,Gamma,G,Vh)...
    + regularize(X,Gamma,'[1/r]',Vh));

DL = 1i/4*(integral(X,Gamma,dG,ntimes(Vh))...
    + regularize(X,Gamma,'grady[1/r]',ntimes(Vh)));