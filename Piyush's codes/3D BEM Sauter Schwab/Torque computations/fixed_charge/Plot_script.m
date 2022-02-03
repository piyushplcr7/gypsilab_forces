loglog(Nvals,abs(torque_mst(:,1)))
hold on
loglog(Nvals,abs(torque_mst(:,2)))
loglog(Nvals,abs(torque_mst(:,3)))

figure
loglog(Nvals,abs(force_mst(:,1)))
hold on
loglog(Nvals,abs(force_mst(:,2)))
loglog(Nvals,abs(force_mst(:,3)))