% Function to find distance of a point from a plane defined through 3
% points
function dist = dist_plane(X,A,B,C)
 
    AB = B-A;
    AC = C-A;

    Normal = cross(AB,AC,2);

%     a = Normal(:,1)
%     b = Normal(:,2)
%     c = Normal(:,3)
    a = Normal(1);
    b = Normal(2);
    c = Normal(3);

    % The plane is ax+by+cz+d=0

%     d = -a.*A(:,1) - b.* A(:,2) - c.* A(:,3)
    d = -a*A(1) - b* A(2) - c* A(3);

%     dist = abs(a.*X(:,1) + b.* X(:,2) + c.* X(:,3) + d)./vecnorm(a.^2+b.^2+c.^2,2,2);

    dist = abs(a*X(1) + b* X(2) + c* X(3) + d)/ sqrt(a^2+b^2+c^2);

end