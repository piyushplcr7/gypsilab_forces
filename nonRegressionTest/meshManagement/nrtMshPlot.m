clean;

%% Point mesh

N = 8;
vtx = randn(N,3); elt = (1:N)'; col = 'b';
m = msh(vtx,elt,col);
figure;
subplot(1,2,1);
plot(m);

axis equal
subplot(1,2,2);
m.col = randn(N,1);
plot(m);
axis equal

%% Segment mesh

N = 30;
m = meshCurve(spirale,N);
m.col = (1 : N)';
figure;
subplot(1,2,1);
plot(m);
axis equal
subplot(1,2,2);
m.col = 'r';
plot(m)
axis equal

%% Triangle mesh

N = 50; rad = 1;
m = mshDisk(N,1);
m.col = 'b';
figure;
subplot(1,2,1)
plot(m);
axis equal
m.col = rand(m.nelt,1);
subplot(1,2,2);
plot(m);
axis equal;

%% Tetrahedral mesh

N = 100; 
m = mshCube(N,[1 1 1]);
figure;
plot(m);
