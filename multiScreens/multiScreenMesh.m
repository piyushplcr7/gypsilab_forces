close all;

ms = nFoldMultiScreen(600,3);
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


Gxy = @(X,Y) femGreenKernel(X,Y,'[1/r]',0);
N = 0;
[~,~,restrictions] = dofMS(phi);

for i = 1:ms.npanels
    Gammai = dom(ms.panels{i},3);
    Ri = restrictions{i};
    phii = fem(ms.panels{i},'P1');
    for j= 1:ms.npanels
        Rj = restrictions{j};
        phij = fem(ms.panels{j},'P1');
        Gammaj = dom(ms.panels{j},3);
        N = N + 1/(4*pi)*Ri'*integral(Gammai,Gammaj,nxgrad(phii),Gxy,nxgrad(phij))*Rj;
        N = N + 1/(4*pi)*Ri'*regularize(Gammai,Gammaj,nxgrad(phii),'[1/r]',nxgrad(phij))*Rj;
    end
end

[P,D] = eig(full(M\N));

[d,ind] = sort(diag(D));
i = find(d> 4,1,'first');
figure;
surf(phi,P(:,ind(i)));
axis equal;
