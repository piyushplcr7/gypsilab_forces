% Function to integrate on the 4 dimensional cube, hard coded n = 5
clear; clc;
f = @(x) sin(x(:,1)+x(:,2)+x(:,3)+x(:,4));
w = [0.568888888888889;
           0.478628670499367;
           0.478628670499367;
           0.236926885056189;
           0.236926885056189];
    x = [0;
            -0.538469310105683;
            0.538469310105683;
            -0.906179845938664;
            0.906179845938664];
    % scaling for integral from 0 to 1
    w = w * 0.5;
    x = (x+1)/2;
    W=[];
    X=[];
    % Tensorize x and w
    %W = kron(kron(w,w'),kron(w,w'));
    %X = kron(kron(x,x'),kron(x,x'));
   
    for i = 1:5
        for j = 1:5
            for k = 1:5
                for l = 1:5
                    W = [W; w(i)*w(j)*w(k)*w(l)];
                    X = [X;x(i), x(j), x(k), x(l)];
                end
            end
        end
    end
    integral = dot(W,f(X))
    %dot(W,f(X))