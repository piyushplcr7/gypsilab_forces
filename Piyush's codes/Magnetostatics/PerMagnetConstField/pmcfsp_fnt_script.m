% Forces and torques script
gpuDevice(1);

delete(gcp('nocreate'));
addpath(genpath("../../../"));

PMCFSP_forcesNtorques(@getMeshSphere,5:11);
PMCFSP_forcesNtorques(@getMeshCube,5:11);
PMCFSP_forcesNtorques(@getMeshCuboid5,5:11);
PMCFSP_forcesNtorques(@getMeshTetraNew,7:11);