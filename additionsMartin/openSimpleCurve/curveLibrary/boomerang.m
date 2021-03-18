function[c] = boomerang()

x = @(t)((0.7 + cos(t*pi)).*cos(t*pi));
y = @(t)(sin(t*pi));
c = SimpleCurve(x,y,[-1,1],'left');

dx = @(t)(-pi*sin(t*pi).*cos(t*pi) - pi*(0.7 + cos(t*pi)).*sin(t*pi));
dy = @(t)(pi*cos(t*pi));
c = supplyDer(c,dx,dy);


end