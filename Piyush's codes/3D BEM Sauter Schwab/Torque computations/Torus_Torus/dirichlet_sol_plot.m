Nvals = Nvals';
sz = size(Nvals,1);

figure;
loglog(Nvals(1:sz-1),abs(averrs(1:sz-1)-averrs(sz)));
title('Av errors');

figure;
loglog(Nvals(1:sz-1),abs(l2errs(1:sz-1)-l2errs(sz)));
title('L2 errors');

figure;
loglog(Nvals(1:sz-1),abs(linferrs(1:sz-1)-linferrs(sz)));
title('Linf errors');

% Checking the pt wise convergence of the trace
figure;
loglog(Nvals(1:sz-1),abs(Tn_plane(1:sz-1)-Tn_plane(sz)));
hold on;
loglog(Nvals(1:sz-1),abs(Tn_nearest(1:sz-1)-Tn_nearest(sz)));
legend('Plane','Nearest DOF');
title('Pt convergence');