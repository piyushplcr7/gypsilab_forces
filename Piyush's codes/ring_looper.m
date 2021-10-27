fileID = fopen('results_ring.txt','a');
fprintf(fileID,'#Fx_bnd Fy_bnd Fx_vol Fy_vol \n',N,F_bnd(1),F_bnd(2),f_vol{1},f_vol{2});
fclose(fileID);
for N = 10:20:250
    forces_ring(N);
end