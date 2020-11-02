function [arc ] = Vcurve


x = @(s)(abs(s));
y = @(s)(s/3);
I = [-1,1];
arc = SimpleCurve(x,y,I);

dx = @(s)(sign(s));
dy = @(s)(1/3 + 0*s);
arc = supplyDer(arc,dx,dy);

end

