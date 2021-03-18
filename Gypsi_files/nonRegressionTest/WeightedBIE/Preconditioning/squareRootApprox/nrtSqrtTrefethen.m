close all;
clear all;
clc;

A = randn(10);
B = A*A';
u = rand(10,1);

C = sqrtm(B);
D = TrefethenSqrt(B,10,[],0.01,100);
norm(C - D)/norm(C)

Du = TrefethenSqrt(B,10,u,0.01,100);
norm(C*u - Du)/norm(C*u)