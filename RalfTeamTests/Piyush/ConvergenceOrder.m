clear all;
close all;

ps = 2:10;

Feng = 0*ps;
Fvol = 0*ps;
FBEM = 0*ps;
errEng = 0*ps;
errVol = 0*ps;
errBEM = 0*ps;
Ns = 0*ps;
Nvol = 0*ps;
Nbound = 0*ps;

refVal = 0.348731673835553;
% refVal = computeNetForce2D(2^(ps(end) +2));

for i = 1:length(ps)
    N = 2^ps(i);
    [FBEM(i),Fvol(i),Feng(i),Nvol(i),Nbound(i),h(i)] = computeNetForce2D(N);
    errEng(i) = abs(Feng(i) - refVal);
    errVol(i) = abs(Fvol(i) - refVal);
    errBEM(i) = abs(FBEM(i) - refVal);
end

loglog(h,errEng);
hold on
loglog(h,errVol);
loglog(h,errBEM);
loglog(h,h*errEng(1)/h(1),'b--');
loglog(h,h.^2*errVol(1)/h(1)^2,'r--');
legend({'Surface formula','Volume formula','BEM formula',...
    'O(h)','O(h^2)'})
xlabel('Mesh size')
ylabel('|Approx - Ref|')
