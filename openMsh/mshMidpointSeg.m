function [mesh,Ir] = mshMidpointSeg(mesh,I)
% Refine segment mesh using midpoint. 

if ~exist('I','var') || isempty(I)
    I = 1:length(mesh);
end

colSave = mesh.col;
mesh = color(mesh,1:length(mesh));


% submesh (segment)
m1 = mesh.sub(I);
m2 = setdiff(mesh,m1);
[A,B] = ABCD(m1);
N = length(m1);
mids = m1.ctr;
vtx = [A;B;mids];
idA = 1:N; idB = N + (1 : N); idMids = 2*N + (1:N);
elt = [[idA(:) idMids(:)];idMids(:) idB(:)];
col = [m1.col;m1.col];
m = msh(vtx,elt,col);
mesh = union(m,m2);

Ir = mesh.col;
mesh.col = colSave(Ir);




end

