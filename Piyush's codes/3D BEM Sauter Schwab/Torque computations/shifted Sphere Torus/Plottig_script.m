data = traces_plane;%vecnorm(tor_torqs,2,2)
N = size(data,1);
figure;
semilogy(Nvals(1:N-1),abs(data(1:N-1) - data(N)));
figure;
loglog(Nvals(1:N-1),abs(data(1:N-1) - data(N)));