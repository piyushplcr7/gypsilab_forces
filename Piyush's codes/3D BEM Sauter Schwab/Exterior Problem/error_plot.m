% Error computation andd plotting script

N = 20:20:160;
av_errs = 0*N;
l2_errs = 0*N;

for i = 1:size(N,2)
    %[l2_errs(i),av_errs(i)] = sphere_torus_exterior_laplace(N(i)); 
    [l2_errs(i),av_errs(i)] = double_sphere_exterior_laplace(N(i));
end

figure;

loglog(N,l2_errs,'o--');
hold on;
loglog(N,av_errs,'o--');
legend('L2 errors','Av errors');