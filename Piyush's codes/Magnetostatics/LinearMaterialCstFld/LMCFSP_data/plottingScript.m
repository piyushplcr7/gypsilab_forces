clear; clc;

addpath(genpath("../../../../"));

prefix = "LMCFSP_forcesNtorques_";

cases = ["getMeshSphere","getMeshCube","getMeshCuboid5","getMeshTetraNew"];


for i = 1:size(cases,2)
        fname = prefix + cases(i)+ ".mat";
        disp(fname);
        load(fname);

        err_Fbem = vecnorm(forces_bem(1:end-1,:)-forces_bem(end,:),2,2);
        err_Fmst = vecnorm(forces_mst(1:end-1,:)-forces_bem(end,:),2,2);
        err_Tbem = vecnorm(torques_bem(1:end-1,:)-torques_bem(end,:),2,2);
        err_Tmst = vecnorm(torques_mst(1:end-1,:)-torques_bem(end,:),2,2);
        h = hvals;

        figure;
        loglog(h(1:end-1),err_Fbem,"+-", 'LineWidth', 2, 'MarkerSize', 10);
        hold on;
        loglog(h(1:end-1),err_Fmst,"^-", 'LineWidth', 2, 'MarkerSize', 10);
        ylabel('Error in total force');
        xlabel('meshwidth h');
        disp(strcat("Fbem rate: ", string(linfitandplot(h(1:end-1),err_Fbem))));
        disp(strcat("Fmst rate: ", string(linfitandplot(h(1:end-1),err_Fmst))));
        legend(["BEM", "MST"],'Location','southeast');
        forces_plot_name = prefix + cases(i)+ "_forces.eps";
        print(forces_plot_name, '-depsc', '-r300');

        figure;
        loglog(h(1:end-1),err_Tbem,"+-", 'LineWidth', 2, 'MarkerSize', 10);
        hold on;
        loglog(h(1:end-1),err_Tmst,"^-", 'LineWidth', 2, 'MarkerSize', 10);
        ylabel('Error in total torque');
        xlabel('meshwidth h');
        disp(strcat("Tbem rate: ", string(linfitandplot(h(1:end-1),err_Tbem))));
        disp(strcat("Tmst rate: ", string(linfitandplot(h(1:end-1),err_Tmst))));
        legend(["BEM", "MST"],'Location','southeast');
        torques_plot_name = prefix + cases(i)+ "_torques.eps";
        print(torques_plot_name, '-depsc', '-r300');
end