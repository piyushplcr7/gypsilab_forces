% Script for TP SP for different meshes
delete(gcp('nocreate'));
addpath(genpath("../../../../"));

gpuDevice(1);
valsrest = 5:12;
valstetra = 7:11;

% TP_VP_dualnorm(@getMeshCuboid5,valsrest);
TP_VP_dualnorm(@getMeshCube,valsrest);
TP_VP_dualnorm(@getMeshSphere,valsrest);
TP_VP_dualnorm(@getMeshTetra1,valstetra);