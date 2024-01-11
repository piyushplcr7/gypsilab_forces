% dualnorm script
gpuDevice(1);

LMCFSP_dualnorm(@getMeshSphere,5:11);
LMCFSP_dualnorm(@getMeshCube,5:11);
LMCFSP_dualnorm(@getMeshCuboid5,5:11);
LMCFSP_dualnorm(@getMeshTetraNew,7:11);