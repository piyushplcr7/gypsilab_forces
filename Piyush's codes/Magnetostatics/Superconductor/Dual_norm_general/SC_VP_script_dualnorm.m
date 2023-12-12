% Script for TP SP for different meshes
delete(gcp('nocreate'));
addpath(genpath("../../../../"));

gpuDevice(2);
valsrest = 5:12;
valstetra = 7:11;

superconductor_VP_dualnorm(@getMeshCuboid5,valsrest);
% superconductor_VP_dualnorm(@getMeshCube,valsrest);
% superconductor_VP_dualnorm(@getMeshSphere,valsrest);
% superconductor_VP_dualnorm(@getMeshTetra1,valstetra);
