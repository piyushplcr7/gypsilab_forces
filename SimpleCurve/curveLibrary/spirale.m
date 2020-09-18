function [arc] = spirale

a=  1;
b= 0.2;
c= 5;
s0 = -0.5;
x = @(s)(a*exp(b*c*(s+s0)).*cos(c*(s+s0)));
y = @(s)(a*exp(b*c*(s+s0)).*sin(c*(s+s0)));
I = [-1,1];
arc = SimpleCurve(x,y,I);

dx = @(s)(a*c*exp(b*c*(s0 + s)).*(b*cos(c*(s0 + s)) - sin(c*(s0 + s))));
dy = @(s)(a*c*exp(b*c*(s0 + s)).*(b*sin(c*(s0 + s)) + cos(c*(s0 + s))));

arc = supplyDer(arc,dx,dy);




end



