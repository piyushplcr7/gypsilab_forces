function out = grad_xyz(X)
    out = [X(:,2).*X(:,3) X(:,1).*X(:,3) X(:,1).*X(:,2)];
end