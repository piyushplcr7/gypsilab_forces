loglog(Nvals,abs(torques_mst(:,1)))
hold on
loglog(Nvals,abs(torques_mst(:,2)))
loglog(Nvals,abs(torques_mst(:,3)))
title('Torques  vs Nvals');

figure
loglog(Nvals,abs(forces_mst(:,1)))
hold on
loglog(Nvals,abs(forces_mst(:,2)))
loglog(Nvals,abs(forces_mst(:,3)))
title('Forces  vs Nvals');

figure
loglog(hvals,abs(torques_mst(:,1)))
hold on
loglog(hvals,abs(torques_mst(:,2)))
loglog(hvals,abs(torques_mst(:,3)))
title('Torques  vs hvals');

figure
loglog(hvals,abs(forces_mst(:,1)))
hold on
loglog(hvals,abs(forces_mst(:,2)))
loglog(hvals,abs(forces_mst(:,3)))
title('Forces  vs hhvals');
