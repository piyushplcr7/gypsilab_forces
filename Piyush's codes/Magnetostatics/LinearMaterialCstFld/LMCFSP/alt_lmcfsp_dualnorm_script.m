% dualnorm script
gpuDevice(1);

delete(gcp('nocreate'));
addpath(genpath("../../../../"));

% LMCFSP_ALT_dualnorm(@getMeshSphere,5:11);
% LMCFSP_ALT_dualnorm(@getMeshCube,5:11);
% LMCFSP_ALT_dualnorm(@getMeshCuboid5,5:11);
% LMCFSP_ALT_dualnorm(@getMeshTetraNew,7:11);

LMCFSP_ALT_dualnormSphere(5:11);

LMCFSP_ALT_dualnormCube(5:11);

LMCFSP_ALT_dualnormCuboid(5:11);