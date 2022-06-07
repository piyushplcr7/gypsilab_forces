% Creating plots
load('tor_tor_combined.mat');

ent = 28;

figure;
fbem = vecnorm(forces_bem,2,2);
loglog(hvals(1:ent-1),abs(fbem(1:ent-1)-fbem(ent))/fbem(ent),'-*','Color','green');
hold on;

fmst = vecnorm(forces_mst,2,2);
loglog(hvals(1:ent-1),abs(fmst(1:ent-1)-fbem(ent))/fbem(ent),'-*','Color','red');

tbem = vecnorm(torques_bem,2,2);
loglog(hvals(1:ent-1),abs(tbem(1:ent-1)-tbem(ent))/tbem(ent),'-+','Color','green');

tmst = vecnorm(torques_mst,2,2);
loglog(hvals(1:ent-1),abs(tmst(1:ent-1)-tbem(ent))/tbem(ent),'-+','Color','red');

legend('Fbem','Fmst','Tbem','Tmst');