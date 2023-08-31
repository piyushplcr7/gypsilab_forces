% Script for PM SP for different meshes
delete(gcp('nocreate'));
addpath(genpath("../../../../"));

gpuDevice(1);

PM_VP_torque(@getMeshCuboid5, 5 : 12, @divcurlfreeM1);
PM_VP_torque(@getMeshCuboid5, 5 : 12, @divcurlfreeM2);

PM_VP_torque(@getMeshSphere, 5:11, @divcurlfreeM1);
PM_VP_torque(@getMeshSphere, 5:11, @divcurlfreeM2);