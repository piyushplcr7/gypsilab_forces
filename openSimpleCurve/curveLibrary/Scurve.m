function [arc ] = Scurve()

a= 0.2;
b= -0.2;
c= 6;
s0 = 0;
x = @(s)(s);
y = @(s)(a*exp(b*c*(s+s0)).*sin(c*(s+s0)));
I = [-1,1];
arc = SimpleCurve(x,y,I);

dx = @(s)(0*s + 1);
dy = @(s)(b*c*a*exp(b*c*(s+s0)).*sin(c*(s+s0)) + c*a*exp(b*c*(s+s0)).*cos(c*(s+s0)));

arc = supplyDer(arc,dx,dy);

end

