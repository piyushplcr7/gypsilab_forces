% Forces and torques script
gpuDevice(1);

LMCFSP_forcesNTorques(@getMeshSphere,5:11);
LMCFSP_forcesNTorques(@getMeshCube,5:11);
LMCFSP_forcesNTorques(@getMeshCuboid5,5:11);
LMCFSP_forcesNTorques(@getMeshTetraNew,7:11);