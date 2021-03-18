%% meshCurve.m
%
% m = meshCurve(c,N) creates a mesh m with N vertices lying on the
% curve c. The nodes are given by x(a + i(b-a)/(N-1)), y(a + i(b-a)/(N-1)) 
% where i = 0..N and x,y is the parametrization of c. 

% Example:

c = Scurve;
N = 50; 
m = meshCurve(c,N);
figure;
plot(m);
% showFrenet(m);

%%
% Using the normal parametrization of the curve, one can obtain a uniform
% mesh with respect to the arclength:

c = spirale; N = 20;
figure 
subplot(1,2,1);
m1 = meshCurve(c,N);
plot(m1);
title('Uniform mesh for intial parametrization')

subplot(1,2,2);
C = normalParam(c);
m2 = meshCurve(C,N);
plot(m2);
title('Uniform mesh for arclength')

%%
% The rep and bds arguments allow to control the density of points. One can
% for example use the curvature as a density parameter for the mesh. 

c = boomerang;
C = normalParam(c);
d = @(t)(sqrt(1 + (2*abs(C.curvature(t)))));
I_d = C.I;

figure
m = meshCurve(C,50,'density',{d,I_d});
plot(m);
title('Mesh-density controlled by curvature')

%%
% One can also choose a change of variables to control the density of the
% mesh. The mesh is uniform with respect to the new parametrization induced
% by this change of variables. 


c = openline(-1,1);
z = @cos; I_z = [-pi,0];
m = meshCurve(c,50,'varChange',{z,I_z});
figure;
plot(m);
title('Chebyshev mesh of (-1,1)')

