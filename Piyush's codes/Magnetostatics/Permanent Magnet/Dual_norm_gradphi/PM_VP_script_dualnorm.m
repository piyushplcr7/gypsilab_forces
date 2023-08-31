% Script for PM SP for different meshes
delete(gcp('nocreate'));
addpath(genpath("../../../../"));

gpuDevice(1);
valsrest = 10 : 11;
valstetra = 7 : 11;

PM_VP_dualnorm_gradphi(@getMeshCuboid5, 10 : 12, @divcurlfreeM1);
PM_VP_dualnorm_gradphi(@getMeshCuboid5, 10 : 12, @divcurlfreeM2);

% PM_VP_dualnorm_gradphi(@getMeshCuboid5, valsrest, @grad_sinx_siny_sinz);
% PM_VP_dualnorm_gradphi(@getMeshCuboid5, valsrest, @grad_xyz);
% PM_VP_dualnorm_gradphi(@getMeshCuboid5, valsrest, @rotConstM);
% PM_VP_dualnorm_gradphi(@getMeshCube, valsrest);
% PM_VP_dualnorm_gradphi(@getMeshSphere, valsrest, @grad_sinx_siny_sinz);

PM_VP_dualnorm_gradphi(@getMeshSphere, valsrest, @divcurlfreeM1);
PM_VP_dualnorm_gradphi(@getMeshSphere, valsrest, @divcurlfreeM2);

% PM_VP_dualnorm_gradphi(@getMeshSphere, valsrest, @grad_xyz);
% PM_VP_dualnorm_gradphi(@getMeshSphere, valsrest, @rotConstM);
% PM_VP_dualnorm_gradphi(@getMeshTetra1, valstetra);