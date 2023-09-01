% Comparing force values from PMSP and PMVP for divcurlfreeM
clear;clc;

load("PMSP_torque_getMeshSphere_divcurlfreeM1.mat");

SPforceBEM = shape_derivatives_bem;
SPforceVOL = shape_derivatives_mst;


load("PMVP_torque_getMeshSphere_divcurlfreeM1.mat");

VPforceBEM = shape_derivatives_bem;
VPforceVOL = shape_derivatives_mst;

disp("M1: SPBEM, SPVOL, VPBEM, VPVOL");
[vecnorm(SPforceBEM,2,2) vecnorm(SPforceVOL,2,2) vecnorm(VPforceBEM,2,2) vecnorm(VPforceVOL,2,2)]

clear;

load("PMSP_torque_getMeshSphere_divcurlfreeM2.mat");

SPforceBEM = shape_derivatives_bem;
SPforceVOL = shape_derivatives_mst;


load("PMVP_torque_getMeshSphere_divcurlfreeM2.mat");

VPforceBEM = shape_derivatives_bem;
VPforceVOL = shape_derivatives_mst;

disp("M2: SPBEM, SPVOL, VPBEM, VPVOL");
[vecnorm(SPforceBEM,2,2) vecnorm(SPforceVOL,2,2) vecnorm(VPforceBEM,2,2) vecnorm(VPforceVOL,2,2)]


