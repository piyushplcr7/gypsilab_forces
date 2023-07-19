function [Vel,DVel] = getSinVelDVel(a,b,c,k)
    
    Vel = @(X) (ones(size(X,1),1)* get_ek(k)).*(sin(a*X(:,1)).* sin(b*X(:,2)).* sin(c*X(:,3)));

    % Get the DVel matrix based on the axis a
    DVel = cell(3,1);

    DVel{1} = @(X) ones(size(X,1),1)* [0 0 0];
    DVel{2} = @(X) ones(size(X,1),1)* [0 0 0];
    DVel{3} = @(X) ones(size(X,1),1)* [0 0 0];

    DVel{k} = @(X) [a*cos(a*X(:,1)).* sin(b*X(:,2)).* sin(c*X(:,3))...
                    b*sin(a*X(:,1)).* cos(b*X(:,2)).* sin(c*X(:,3))...
                    c*sin(a*X(:,1)).* sin(b*X(:,2)).* cos(c*X(:,3))];
end