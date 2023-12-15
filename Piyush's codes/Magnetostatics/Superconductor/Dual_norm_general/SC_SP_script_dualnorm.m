% Script for TP SP for different meshes
delete(gcp('nocreate'));
addpath(genpath("../../../../"));

gpuDevice(3);
valsrest = 5:12;
valstetra = 7:11;

% superconductor_SP_dualnorm(@getMeshCuboid5,valsrest);
% superconductor_SP_dualnorm(@getMeshCube,valsrest);
% superconductor_SP_dualnorm(@getMeshSphere,valsrest);
% superconductor_SP_dualnorm(@getMeshTetra1,valstetra);
superconductor_SP_dualnorm(@getMeshTetraNew,valstetra);
