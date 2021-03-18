function [ curve ] = circle(c,r)
% creates the object of class Curve corresponding to the circle C(r,c),
% where c is the center and r is the radius 

if nargin == 0
    c = [0,0];
    r = 1;  
end

x = @(t)(c(1) + r*cos(t));
y = @(t)(c(2) + r*sin(t));
I = [-pi,pi];
boundedSide = 'left';
curve = SimpleCurve(x,y,I,boundedSide);

curve = supplyDer(curve,@(t)(-r*sin(t)),@(t)(r*cos(t)));

end

