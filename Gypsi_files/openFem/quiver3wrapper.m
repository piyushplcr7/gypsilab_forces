function a = quiver3wrapper(X,Y,col)
    quiver3(X(:,1),X(:,2),X(:,3),Y(:,1),Y(:,2),Y(:,3),col);
    a = [];
    return;
end