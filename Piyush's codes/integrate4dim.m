% Function to integrate on the 4 dimensional cube, hard coded n = 5
function integral = integrate4dim(f)
%load('X','X');
%load('W','W');
global X;
global W;
%     w = [0.568888888888889;
%            0.478628670499367;
%            0.478628670499367;
%            0.236926885056189;
%            0.236926885056189];
%     x = [0;
%             -0.538469310105683;
%             0.538469310105683;
%             -0.906179845938664;
%             0.906179845938664];
%     % scaling for integral from 0 to 1
%     w = w * 0.5;
%     x = (x+1)/2;
    
%     integral2 = 0;
    integral = dot(W,f(X(:,1),X(:,2),X(:,3),X(:,4)));
%    for i = 1:5
%        for j = 1:5
%            for k = 1:5
%                for l = 1:5
%                    integral2 = integral2 + w(i)*w(j)*w(k)*w(l)*f(x(i),x(j),x(k),x(l));
%                end
%            end
%        end
%    end
   % disp(integral-integral2);
end