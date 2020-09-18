function[arc] = semicircle(c,r)

if nargin == 0
    c = [0,0];
    if nargin <= 1
        r = 1;
    end
end

arc = circle(c,r);
arc = portion(arc,[0,pi]);



end