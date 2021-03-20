close all;

ms = nFoldMultiScreen(500,3);
figure 
plot(ms);
hold on
plotNrm(ms.panels{1},'r');
plotNrm(ms.panels{2},'g');
plotNrm(ms.panels{3},'b');
plot(ms.bnd,'y');

phi = msFem(ms,'multi_P1');
Gamma = msDom(ms,3);

M = integral(Gamma,phi,phi);
K = integral(Gamma,grad(phi),grad(phi));
[P,D] = eig(full(M\K));

[d,ind] = sort(diag(D));
figure;
surf(phi,P(:,ind(4)));

axis equal