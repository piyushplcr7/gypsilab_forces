% Script for PM SP for different meshes
delete(gcp('nocreate'));
addpath(genpath("../../../../"));

gpuDevice(1);
valsrest = 5:12;
valstetra = 7:11;

PM_VP_dualnorm_mst(@getMeshCuboid5,valsrest);
PM_VP_dualnorm_mst(@getMeshCube,valsrest);
PM_VP_dualnorm_mst(@getMeshSphere,valsrest);
PM_VP_dualnorm_mst(@getMeshTetra1,valstetra);