% Comparing force values from PMSP and PMVP for divcurlfreeM
clear;clc;

load("PMSP_dualnorm_getMeshSphere_divcurlfreeM1.mat");

SPforceBEM = [shape_derivatives_bem(:,1) shape_derivatives_bem(:,28) shape_derivatives_bem(:,55)];
SPforceVOL = [shape_derivatives_mst(:,1) shape_derivatives_mst(:,28) shape_derivatives_mst(:,55)];


load("PMVP_dualnorm_getMeshSphere_divcurlfreeM1.mat");

VPforceBEM = [shape_derivatives_bem(:,1) shape_derivatives_bem(:,28) shape_derivatives_bem(:,55)];
VPforceVOL = [shape_derivatives_mst(:,1) shape_derivatives_mst(:,28) shape_derivatives_mst(:,55)];

disp("M1: SPBEM, SPVOL, VPBEM, VPVOL");
[vecnorm(SPforceBEM,2,2) vecnorm(SPforceVOL,2,2) vecnorm(VPforceBEM,2,2) vecnorm(VPforceVOL,2,2)]

clear;

load("PMSP_dualnorm_getMeshSphere_divcurlfreeM2.mat");

SPforceBEM = [shape_derivatives_bem(:,1) shape_derivatives_bem(:,28) shape_derivatives_bem(:,55)];
SPforceVOL = [shape_derivatives_mst(:,1) shape_derivatives_mst(:,28) shape_derivatives_mst(:,55)];


load("PMVP_dualnorm_getMeshSphere_divcurlfreeM2.mat");

VPforceBEM = [shape_derivatives_bem(:,1) shape_derivatives_bem(:,28) shape_derivatives_bem(:,55)];
VPforceVOL = [shape_derivatives_mst(:,1) shape_derivatives_mst(:,28) shape_derivatives_mst(:,55)];

disp("M2: SPBEM, SPVOL, VPBEM, VPVOL");
[vecnorm(SPforceBEM,2,2) vecnorm(SPforceVOL,2,2) vecnorm(VPforceBEM,2,2) vecnorm(VPforceVOL,2,2)]


