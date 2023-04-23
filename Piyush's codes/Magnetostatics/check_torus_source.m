% Check torus surface source
clear;clc;

[J,mesh] = get_torus_source(500,3,0.5);
omega = dom(mesh,3);

[X,W] = omega.qud;

JX = J(X(:,1),X(:,2),X(:,3));

plot(mesh);
hold on;
quiver3(X(:,1),X(:,2),X(:,3),JX(:,1),JX(:,2),JX(:,3));