clear; clc;

addpath(genpath("../../../../"));

processDualNormData('LMCFVP_dualnorm_getMeshCube.mat',1:7,1:6,1:7);

processDualNormData('LMCFVP_dualnorm_getMeshSphere.mat',1:7,1:6,1:7);

processDualNormData('LMCFVP_dualnorm_getMeshCuboid5.mat',1:7,1:6,1:7);

processDualNormData('LMCFVP_dualnorm_getMeshTetraNew.mat',1:5,1:4,1:5);