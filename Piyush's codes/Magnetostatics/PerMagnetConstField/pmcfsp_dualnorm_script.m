% dualnorm script
gpuDevice(1);

delete(gcp('nocreate'));
addpath(genpath("../../../"));

PMCFSP_dualnorm(@getMeshSphere,5:11);
PMCFSP_dualnorm(@getMeshCube,5:11);
PMCFSP_dualnorm(@getMeshCuboid5,5:11);
PMCFSP_dualnorm(@getMeshTetraNew,7:11);