% dualnorm script
gpuDevice(1);

DielSC_dualnorm(@getMeshSphere,5:11);
DielSC_dualnorm(@getMeshCube,5:11);
DielSC_dualnorm(@getMeshCuboid5,5:11);
DielSC_dualnorm(@getMeshTetraNew,7:11);