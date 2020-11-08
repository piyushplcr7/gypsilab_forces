clc;

disp(char([087 097 115 032 108 097 117 102 116 032 117 110 100 032 119 097 115 032 119 105 114 100 032 103 101 119 111 114 102 101 110 063 013 010]));
addpath(genpath(pwd));

run nonRegressionTest/WeightedBIE/Preconditioning/nrtPrecHelmholtzCurves.m;