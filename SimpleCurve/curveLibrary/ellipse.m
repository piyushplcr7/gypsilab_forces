function[c] = ellipse(a,b)


x = @(t)(a*cos(pi*t));
y = @(t)(b*sin(pi*t));
I = [-1 1];

c = SimpleCurve(x,y,I,'left');

dx = @(t)(-pi*a*sin(pi*t));
dy = @(t)(pi*b*cos(pi*t));

c = supplyDer(c,dx,dy);

end