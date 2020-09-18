function[c] = parabola

x = @(t)(t);
y = @(t)(t.^2 - 0.25);
I = [-1,1];
dx = @(t)(0*t + 1); dy = @(t)(2*t);
c = supplyDer(SimpleCurve(x,y,I),dx,dy);


end