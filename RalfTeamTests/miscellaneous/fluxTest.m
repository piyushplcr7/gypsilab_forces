m = mshCube(20,[1 1 1]);
Vh = fem(m,'RWG');
uh = ones(Vh.ndof,1);
traceU = normalTrace_RT(Vh)*uh;
plot(dm,traceU)
hold on
plotNrm(dm,'r')
pltNaturalNrm(dm,'g');