% Forces and torques script
gpuDevice(1);

delete(gcp('nocreate'));
addpath(genpath("../../../"));

DielSC_forcesNtorques(@getMeshSphere,5:11);
DielSC_forcesNtorques(@getMeshCube,5:11);
DielSC_forcesNtorques(@getMeshCuboid5,5:11);
DielSC_forcesNtorques(@getMeshTetraNew,7:11);