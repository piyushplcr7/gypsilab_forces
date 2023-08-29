function out = grad_sinx_siny_sinz(X)
    out = [cos(X(:,1)).*sin(X(:,2)).*sin(X(:,3))...
           sin(X(:,1)).*cos(X(:,2)).*sin(X(:,3))...
           sin(X(:,1)).*sin(X(:,2)).*cos(X(:,3))];
end