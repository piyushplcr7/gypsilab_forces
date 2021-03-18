function [arc ] = Vcurve


y = @(s)(abs(s));
x = @(s)(s);
I = [-1,1];
arc = SimpleCurve(x,y,I);

dy = @(s)(sign(s));
dx = @(s)(0*s + 1);
arc = supplyDer(arc,dx,dy);

end

