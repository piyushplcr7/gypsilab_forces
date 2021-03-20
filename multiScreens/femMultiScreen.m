function [P] = femMultiScreen(Xfem,Xdom)
% Dof to quadrature matrix

ms = Xdom.msh;
gs = Xdom.gss;
[~,W,I] = qud(Xdom);
[~,~,restrictions] = dofMS(Xfem);
pani = ms.panels{1};
domi = dom(pani,gs);
femi = fem(pani,Xfem.typ(7:end));
femi.opr = Xfem.opr;
uqmi = femi.uqm(domi);
idx = I{1};
jdx = 1:size(domi.qud,1);
Q = sparse(idx,jdx,1,length(W),length(jdx));
if iscell(uqmi)
    P = cell(length(uqmi),1);
    for d = 1:length(uqmi)
        P{d} = Q*uqmi{d}*restrictions{1};
    end
else
    P =  Q*uqmi*restrictions{1};
end

for i = 2:ms.npanels
    pani = ms.panels{i};
    domi = dom(pani,gs);
    femi = fem(pani,Xfem.typ(7:end));
    femi.opr = Xfem.opr;
    uqmi = femi.uqm(domi);
    idx = I{i};
    jdx = 1:size(domi.qud,1);
    Q = sparse(idx,jdx,1,length(W),length(jdx));
    if iscell(uqmi)
        for d = 1:length(uqmi)
            P{d} = P{d} + Q*uqmi{d}*restrictions{i};
        end
    else
        P =  P + Q*uqmi*restrictions{i};
    end
end
end