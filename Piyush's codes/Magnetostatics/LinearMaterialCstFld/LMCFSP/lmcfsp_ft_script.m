% Forces and torques script
gpuDevice(1);

delete(gcp('nocreate'));
addpath(genpath("../../../../"));

LMCFSP_forcesNtorques(@getMeshSphere,5:11);
LMCFSP_forcesNtorques(@getMeshCube,5:11);
LMCFSP_forcesNtorques(@getMeshCuboid5,5:11);
LMCFSP_forcesNtorques(@getMeshTetraNew,7:11);