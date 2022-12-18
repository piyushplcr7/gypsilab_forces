u = @(x,y,z) [cos(pi*x).*sin(pi*y).*sin(pi*z) -0.5*sin(pi*x).*cos(pi*y).*sin(pi*z) -0.5*sin(pi*x).*sin(pi*y).*cos(pi*z)];

t = 0:0.05:1;

[X,Y,Z] = meshgrid(t,t,t);