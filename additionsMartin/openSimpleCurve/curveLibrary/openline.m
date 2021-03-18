function [ curve ] = openline(a,b)

x = @(t)(a + (t+1)/2*(b-a));
y = @(t)(0*t);
I = [-1,1];
curve = SimpleCurve(x,y,I);

dx= @(t)((b-a)/2 + 0*t);
dy = @(t)(0*t);
curve = supplyDer(curve,dx,dy);


end

