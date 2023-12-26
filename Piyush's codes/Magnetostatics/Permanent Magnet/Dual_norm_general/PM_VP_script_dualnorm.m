% Script for PM SP for different meshes
delete(gcp('nocreate'));
addpath(genpath("../../../../"));

gpuDevice(1);
valsrest = 5:12;
valstetra = 7:11;

% PM_VP_dualnorm(@getMeshCuboid5,valsrest);
% PM_VP_dualnorm(@getMeshCube,valsrest);
% PM_VP_dualnorm(@getMeshSphere,valsrest);
PM_VP_dualnorm(@getMeshTetraNew,valstetra);