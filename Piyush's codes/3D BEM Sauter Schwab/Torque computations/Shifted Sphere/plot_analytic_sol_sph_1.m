Nvals = Nvals';
sz = size(Nvals,1);

figure;
loglog(Nvals,averrs);
title('Av errors');

figure;
loglog(Nvals,l2errs);
title('L2 errors');

figure;
loglog(Nvals,linferrs);
title('Linf errors');

% Checking the pt wise convergence of the trace
figure;
loglog(Nvals(1:sz-1),abs(Tn_plane(1:sz-1)-Tn_plane(sz)));
hold on;
loglog(Nvals(1:sz-1),abs(Tn_nearest(1:sz-1)-Tn_nearest(sz)));
legend('Plane','Nearest DOF');
title('Pt convergence');

%%
% Checking the pt wise convergence of the trace
figure;
loglog(Nvals(1:sz-1),abs(Tn_plane(1:sz-1)+1/Rad));
hold on;
loglog(Nvals(1:sz-1),abs(Tn_nearest(1:sz-1)+1/Rad));
legend('Plane','Nearest DOF');
title('Pt convergence');