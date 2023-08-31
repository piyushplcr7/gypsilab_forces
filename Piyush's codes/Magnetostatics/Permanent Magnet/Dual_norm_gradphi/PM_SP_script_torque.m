% Script for PM SP for different meshes
delete(gcp('nocreate'));
addpath(genpath("../../../../"));

gpuDevice(1);
valsrest = 10 : 11;
valstetra = 7 : 11;


PM_SP_torque(@getMeshCuboid5, 5 : 12, @divcurlfreeM1);
PM_SP_torque(@getMeshCuboid5, 5 : 12, @divcurlfreeM2);

PM_SP_torque(@getMeshSphere, 5:11, @divcurlfreeM1);
PM_SP_torque(@getMeshSphere, 5:11, @divcurlfreeM2);
