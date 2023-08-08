function [Vel,DVel] = getCosVelDVel(a,b,c,alpha)
    
    Vel = @(X) (ones(size(X,1),1)* get_ek(alpha)).*(cos(a*X(:,1)).* cos(b*X(:,2)).* cos(c*X(:,3)));

    % Get the DVel matrix based on the axis a
    DVel = cell(3,1);

    DVel{1} = @(X) ones(size(X,1),1)* [0 0 0];
    DVel{2} = @(X) ones(size(X,1),1)* [0 0 0];
    DVel{3} = @(X) ones(size(X,1),1)* [0 0 0];

    DVel{alpha} = @(X) [-a*sin(a*X(:,1)).* cos(b*X(:,2)).* cos(c*X(:,3))...
                    -b*cos(a*X(:,1)).* sin(b*X(:,2)).* cos(c*X(:,3))...
                    -c*cos(a*X(:,1)).* cos(b*X(:,2)).* sin(c*X(:,3))];
end