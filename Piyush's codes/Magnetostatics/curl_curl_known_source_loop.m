% Curl curl known source loop
clear; clc;
addpath(genpath("../../"));
for N = 2:12
    disp(2^N);
    [l2err,~,~] = eval_err_known_source(2^N)
end