function [Vel,DVel] = getPolyVelDVel(a,b,c,k)
    
    Vel = @(X) (ones(size(X,1),1)* get_ek(k)).*(X(:,1).^a .* X(:,2).^b .* X(:,3).^c);

    % Get the DVel matrix based on the axis a
    DVel = cell(3,1);

    DVel{1} = @(X) ones(size(X,1),1)* [0 0 0];
    DVel{2} = @(X) ones(size(X,1),1)* [0 0 0];
    DVel{3} = @(X) ones(size(X,1),1)* [0 0 0];

    DVel{k} = @(X) [a*X(:,1).^(a-1).*X(:,2).^b.*X(:,3).^c...
                    b*X(:,1).^a .* X(:,2).^(b-1) .* X(:,3).^c...
                    c*X(:,1).^a .* X(:,2).^b .* X(:,3).^(c-1)];
end