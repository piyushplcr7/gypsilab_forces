% Forces and torques script
delete(gcp('nocreate'));
addpath(genpath("../../../"));

LMCFVP_forcesNtorques(@getMeshSphere,5:11);
LMCFVP_forcesNtorques(@getMeshCube,5:11);
LMCFVP_forcesNtorques(@getMeshCuboid5,5:11);
LMCFVP_forcesNtorques(@getMeshTetraNew,7:11);