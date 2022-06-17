mesh = mshDisk(10,1);

V = fem(mesh,'P1');

omega = dom(mesh,3);

uqmat = V.uqm(omega);

[X,W,elt2qud] = omega.qud;

% Visualizing basis function
plot(mesh);
hold on;
scatter3(X(:,1),X(:,2),uqmat(:,2));

