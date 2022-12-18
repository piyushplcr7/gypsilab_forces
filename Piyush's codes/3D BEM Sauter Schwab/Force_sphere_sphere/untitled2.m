clear;clc;
N = 1:9000;

psi = 0*N;

for n=N
    psi(n) = log_derv_gamma(1,n);
end