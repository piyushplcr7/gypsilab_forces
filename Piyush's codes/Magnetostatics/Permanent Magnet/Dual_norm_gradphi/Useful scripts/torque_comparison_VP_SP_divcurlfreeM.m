%% Sphere M1
clear;clc;

load("PMSP_torque_getMeshSphere_divcurlfreeM1.mat");

SPforceBEM = shape_derivatives_bem;
SPforceVOL = shape_derivatives_mst;


load("PMVP_torque_getMeshSphere_divcurlfreeM1.mat");

VPforceBEM = shape_derivatives_bem;
VPforceVOL = shape_derivatives_mst;

disp("M1: SPBEM, SPVOL, VPBEM, VPVOL");
[vecnorm(SPforceBEM,2,2) vecnorm(SPforceVOL,2,2) vecnorm(VPforceBEM,2,2) vecnorm(VPforceVOL,2,2)]

%% Sphere M2
clear;

load("PMSP_torque_getMeshSphere_divcurlfreeM2.mat");

SPforceBEM = shape_derivatives_bem;
SPforceVOL = shape_derivatives_mst;


load("PMVP_torque_getMeshSphere_divcurlfreeM2.mat");

VPforceBEM = shape_derivatives_bem;
VPforceVOL = shape_derivatives_mst;

disp("M2: SPBEM, SPVOL, VPBEM, VPVOL");
[vecnorm(SPforceBEM,2,2) vecnorm(SPforceVOL,2,2) vecnorm(VPforceBEM,2,2) vecnorm(VPforceVOL,2,2)]


%% Cuboid M1

clear;clc;

load("PMSP_torque_getMeshCuboid5_divcurlfreeM1.mat");

SPforceBEM = shape_derivatives_bem;
SPforceVOL = shape_derivatives_mst;


load("PMVP_torque_getMeshCuboid5_divcurlfreeM1.mat");

VPforceBEM = shape_derivatives_bem;
VPforceVOL = shape_derivatives_mst;

disp("M1: SPBEM, SPVOL, VPBEM, VPVOL");
[vecnorm(SPforceBEM,2,2) vecnorm(SPforceVOL,2,2) vecnorm(VPforceBEM,2,2) vecnorm(VPforceVOL,2,2)]

%% Cuboid M2
clear;clc;

load("PMSP_torque_getMeshCuboid5_divcurlfreeM2.mat");

SPforceBEM = shape_derivatives_bem;
SPforceVOL = shape_derivatives_mst;


load("PMVP_torque_getMeshCuboid5_divcurlfreeM2.mat");

VPforceBEM = shape_derivatives_bem;
VPforceVOL = shape_derivatives_mst;

disp("M1: SPBEM, SPVOL, VPBEM, VPVOL");
[vecnorm(SPforceBEM,2,2) vecnorm(SPforceVOL,2,2) vecnorm(VPforceBEM,2,2) vecnorm(VPforceVOL,2,2)]

%% Torque convergence plot

smallsize = min(size(VPforceBEM,1),size(SPforceBEM,1));

if size(VPforceBEM,1) == smallsize
    converged = SPforceBEM(end,:);
else
    converged = VPforceBEM(end,:);
end

loglog(hvals(1:smallsize),vecnorm(SPforceBEM(1:smallsize,:)-converged,2,2));
hold on;
loglog(hvals(1:smallsize),vecnorm(VPforceBEM(1:smallsize,:)-converged,2,2));
loglog(hvals(1:smallsize),vecnorm(SPforceVOL(1:smallsize,:)-converged,2,2));
loglog(hvals(1:smallsize),vecnorm(VPforceVOL(1:smallsize,:)-converged,2,2));

legend(["SPBEM","VPBEM","SPVOL","VPVOL"]);



