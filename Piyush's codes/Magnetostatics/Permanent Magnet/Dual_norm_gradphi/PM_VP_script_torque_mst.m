% Script for PM SP for different meshes
delete(gcp('nocreate'));
addpath(genpath("../../../../"));

% gpuDevice(1);

PM_VP_torque_mst(@getMeshCuboid5, 5 : 12, @divcurlfreeM1);
PM_VP_torque_mst(@getMeshCuboid5, 5 : 12, @divcurlfreeM2);

PM_VP_torque_mst(@getMeshSphere, 5:12, @divcurlfreeM1);
PM_VP_torque_mst(@getMeshSphere, 5:12, @divcurlfreeM2);