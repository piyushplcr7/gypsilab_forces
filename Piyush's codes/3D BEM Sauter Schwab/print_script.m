% Printing script cube cube
clear;clc;close all;
forces_bnd_in = [];
forces_bnd_out = [];
forces_sg = [];
h = [0.055555555555556; 0.055555555555556;  0.031250000000000; 0.020000000000000; 0.013888888888889; 0.007812500000000; 0.005000000000000; 0.002958579881657];

base = "Cube_Cube_";

% Extracting the data
for i=1:8
    N = 2^(i+3)
    fname = append(base,int2str(N));
    
    % Loading the variables
    load(fname,'forcex');
    load(fname,'forcey');
    load(fname,'forcez');

    load(fname,'classical_force_in');
    load(fname,'classical_force_out');
    
    % Saving the variables at the correct location
    forces_bnd_in = [forces_bnd_in;  norm(classical_force_in)];
    forces_bnd_out = [forces_bnd_out;  norm(classical_force_out)];
    forces_sg = [forces_sg; norm([forcex, forcey, forcez])];
end

% Splitting the data into data and converged value
f_conv = forces_sg(8);

fbndiner = abs(forces_bnd_in(1:7)-f_conv);
fbndouter = abs(forces_bnd_out(1:7)-f_conv);
fsger = abs(forces_sg(1:7)-f_conv);

h = h(1:7);

loglog(h,fbndiner,'-or');
hold on;
loglog(h,fbndouter,'-*g');
loglog(h,fsger,'-+b');
legend('Boundary in','Boundary out','BEM','Location','southeast');
xlabel('Meshwidth h');
ylabel('Absolute error');
title('Cube Cube');
xlim([4e-3 7e-2]);
ylim([8e-5 2e-1]);