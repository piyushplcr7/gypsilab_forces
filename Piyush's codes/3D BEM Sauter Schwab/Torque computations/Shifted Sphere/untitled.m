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
exact_val = -1/4/pi/36;
figure;
%loglog(Nvals(1:sz-1),abs(Tn_plane(1:sz-1)-Tn_plane(sz)));
loglog(Nvals,abs(Tn_plane-exact_val));
hold on;
loglog(Nvals(1:sz),abs(Tn_nearest-exact_val));
legend('Plane','Nearest DOF');
title('Pt convergence');