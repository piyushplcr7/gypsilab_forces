% Forces and torques script
gpuDevice(1);
delete(gcp('nocreate'));
addpath(genpath("../../../../"));

% LMCFVP_dualnorm(@getMeshSphere,5:11);
% LMCFVP_dualnorm(@getMeshCube,5:11);
% LMCFVP_dualnorm(@getMeshCuboid5,5:11);
% LMCFVP_dualnorm(@getMeshTetraNew,7:11);

LMCFVP_dualnorm(@sphSphMesh,5:11);
LMCFVP_dualnorm(@cubeSphMesh,5:11);
LMCFVP_dualnorm(@cuboidSphMesh,5:11);
LMCFVP_dualnorm(@tetraSphMesh,7:11);