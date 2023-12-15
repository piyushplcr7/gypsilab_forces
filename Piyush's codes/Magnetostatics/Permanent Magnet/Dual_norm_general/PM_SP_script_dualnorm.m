% Script for PM SP for different meshes
delete(gcp('nocreate'));
addpath(genpath("../../../../"));

gpuDevice(2);
valsrest = 5:12;
valstetra = 7:11;

% PM_SP_dualnorm(@getMeshCuboid5,valsrest);
% PM_SP_dualnorm(@getMeshCube,valsrest);
% PM_SP_dualnorm(@getMeshSphere,valsrest);
PM_SP_dualnorm(@getMeshTetraNew,valstetra);
