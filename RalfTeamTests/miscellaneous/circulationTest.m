m = mshCube(20,[1 1 1]);
Vh = fem(m,'NED');
uh = ones(Vh.ndof,1);
traceU = tangentialTrace_NED(Vh)*uh;
plot(dm,traceU)
hold on
plotNrm(dm,'r')
pltNaturalNrm(dm,'g');