clear all;
close all;
m = mshAssymetricRing(2,1,2);
S1 = fem(m,'P1');
S0 = fem(m,'P0');
S1.ndof
S0.ndof
s0dofs = S0.dof;
s1dofs = S1.dof;

% Visualizing the dofs
plot(m);
hold on;
scatter(s0dofs(:,1),s0dofs(:,2),'red');
scatter(s1dofs(:,1),s1dofs(:,2),'blue');
