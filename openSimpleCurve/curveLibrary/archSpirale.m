function [arc] = archSpirale()

m = 3;
y = @(s)((1 + 0.5*m*pi*(s + 1)).*cos(m*pi*(s+1))/5);
x = @(s)(-(1 + 0.5*m*pi*(s + 1)).*sin(m*pi*(s+1))/5);

I = [-1,1];
arc = SimpleCurve(x,y,I);

dy = @(s)(-m*pi*(1 + 0.5*m*pi*(s + 1)).*sin(m*pi*(s+1))/5 ...
    + 0.5*m*pi.*cos(m*pi*(s+1))/5);
dx =  @(s)(-m*pi*(1 + 0.5*m*pi*(s + 1)).*cos(m*pi*(s+1))/5 ...
    -0.5*m*pi*sin(m*pi*(s+1))/5);

arc = supplyDer(arc,dx,dy);


end

