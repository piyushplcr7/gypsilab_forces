function[c] = parabola

x = @(t)(t);
y = @(t)(1/2*t.^2 - 1);
I = [-1,1];
dx = @(t)(0*t + 1); dy = @(t)(t);
c = supplyDer(SimpleCurve(x,y,I),dx,dy);


end