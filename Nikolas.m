clc;

disp(char([086 105 101 108 101 110 032 068 097 110 107 044 032 115 101 104 114 032 108 105 101 098 101 114 032 078 105 107 111 108 097 115 033]))
addpath(genpath(pwd));

pause(0.8);
run nonRegressionTest/WeightedBIE/Preconditioning/nrtPrecHelmholtzCurves.m;