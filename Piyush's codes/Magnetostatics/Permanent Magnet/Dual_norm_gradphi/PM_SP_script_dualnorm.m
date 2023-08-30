% Script for PM SP for different meshes
delete(gcp('nocreate'));
addpath(genpath("../../../../"));

gpuDevice(1);
valsrest = 5:9;
valstetra = 7:11;

% PM_SP_dualnorm_gradphi(@getMeshCuboid5,valsrest);
% PM_SP_dualnorm_gradphi(@getMeshCube,valsrest);
% PM_SP_dualnorm_gradphi(@getMeshSphere,valsrest);
% PM_SP_dualnorm_gradphi(@getMeshTetra1,valstetra);

PM_SP_dualnorm_gradphi(@getMeshCuboid5,valsrest,@divcurlfreeM1);
PM_SP_dualnorm_gradphi(@getMeshCuboid5,valsrest,@divcurlfreeM2);

PM_SP_dualnorm_gradphi(@getMeshCuboid5,valsrest,@divcurlfreeM1);
PM_SP_dualnorm_gradphi(@getMeshCuboid5,valsrest,@divcurlfreeM2);