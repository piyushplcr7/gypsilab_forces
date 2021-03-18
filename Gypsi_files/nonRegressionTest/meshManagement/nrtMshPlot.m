clean;

%% 1°) PLOT MESH ITSLEF


%% Point mesh

N = 8;
vtx = randn(N,3); elt = (1:N)'; col = 'b';
m = msh(vtx,elt,col);
disp(m);
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
disp(m);
figure;
subplot(1,2,1);
plot(m);
hold on
plotNrm(m);
plotTgt(m);
axis equal
subplot(1,2,2);
m.col = 'r';
m = swap(m);
plot(m)
hold on
plotUnitNrm(m,'g');
hold on
plotUnitTgt(m,'b');
axis equal

%% Triangle mesh

N = 50; rad = 1;
m = mshDisk(N,1);
m.vtx(:,3) = sqrt(1 - norm2D(m.vtx(:,1:2)));
disp(m);
figure;
subplot(1,2,1)
plot(m,'b');
hold on
plotTgt(m)
axis equal
view(178,28);
m.col = rand(m.nelt,1);
subplot(1,2,2);
plot(m);
hold on
plotNrm(m);
axis equal;
view(178,28);

%% Tetrahedral mesh

N = 8; 
m = mshCube(N,[1 1 1]);
m.col = (1:m.nelt)';
disp(m);
figure;
subplot(1,2,1);
plot(m);
axis equal;
view(12,46);
subplot(1,2,2);
plot(m.explode(0.85));
hold on
plotTgt(m.explode(0.85),'r');
plotNrm(m.explode(0.85).fce,'r');
axis equal;
view(12,46);
% We can't plot normals on a tetrahedral mesh.


%% 2°) Plot functions defined on the mesh


%% Point mesh

m = prt(mshSphere(50,1)); 
f = @(X)(cos(X(:,1) + sin(X(:,2))));
figure;
plot(m,f)

%% Segment mesh


m = meshCurve(openline(-1,1),10); 
f = @(X)(cos(X(:,1)));
plotOn(m,f);

%% 


