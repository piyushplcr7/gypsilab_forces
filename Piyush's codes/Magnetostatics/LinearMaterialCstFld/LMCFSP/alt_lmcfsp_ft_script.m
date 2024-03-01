% Forces and torques script

delete(gcp('nocreate'));
addpath(genpath("../../../../"));

%LMCFSP_forcesNtorquesAlt(@getMeshSphere,5:11);
%LMCFSP_forcesNtorquesAlt(@getMeshCube,5:11);
%LMCFSP_forcesNtorquesAlt(@getMeshCuboid5,5:11);
% LMCFSP_forcesNtorquesAlt(@getMeshTetraNew,7:11);

LMCFSP_forcesNtorquesAltSphere(5:11);

LMCFSP_forcesNtorquesAltCube(5:11);

LMCFSP_forcesNtorquesAltCuboid(5:11);