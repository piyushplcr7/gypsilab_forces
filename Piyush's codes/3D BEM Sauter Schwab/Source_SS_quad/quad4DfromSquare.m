function [Xss, Wss] = quad4DfromSquare(rule)

[xsq,wsq] = rule();

% Rescaling from -1,1 to 0,1
wsq = wsq/4;
xsq = (1+xsq)/2;

wsq = wsq';
xsq = xsq';


nq = size(wsq, 1);

X = zeros(nq^2, 4);
W = zeros(nq^2, 1);

count = 1;
for i1 = 1:nq
   for i2 = 1:nq
          X(count, :) = [xsq(i1,:) xsq(i2,:)];
          W(count)    = wsq(i1) * wsq(i2);
          count = count + 1;
    end
end

% Don't know why I need to scale like this
W = W * 2;
   

[Xss{1}, Wss{1}] = quad1_bndry(X, W); % identical

[Xss{2}, Wss{2}] = quad2_bndry(X, W); % common edge

[Xss{3}, Wss{3}] = quad3_bndry(X, W); % common vertex

[Xss{4}, Wss{4}] = quad4_bndry(X, W); % far away


end