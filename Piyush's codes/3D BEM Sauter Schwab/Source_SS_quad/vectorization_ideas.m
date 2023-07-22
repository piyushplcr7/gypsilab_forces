%Elt = sparse(bndmesh.Nelt,bndmesh.nvtx);
elts = 1:bndmesh.nelt;
elts = elts';
elts = repelem(elts,3);
vtcs = reshape(bndmesh.elt',[bndmesh.nelt*3 1]);

Eltmat = sparse(elts,vtcs,ones(size(vtcs)),bndmesh.nelt,bndmesh.nvtx);

Intmat = Eltmat * Eltmat';

% othermat = sparse(bndmesh.nelt,bndmesh.nelt);
% Create another way to verify % Verified!
% for i = 1:bndmesh.nelt
%     for j = 1:bndmesh.nelt
%         intersection = intersect(bndmesh.elt(i,:),bndmesh.elt(j,:));
%         l = length(intersection);
%         othermat(i,j) = l;
% 
%     end
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% FAR AWAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[I0,J0] = find(Intmat == 0);

elts_I0 = bndmesh.elt(I0,:);
elts_J0 = bndmesh.elt(J0,:);

AI0 = bndmesh.vtx(elts_I0(1,:),:);
BI0 = bndmesh.vtx(elts_I0(2,:),:);
CI0 = bndmesh.vtx(elts_I0(3,:),:);

AJ0 = bndmesh.vtx(elts_J0(1,:),:);
BJ0 = bndmesh.vtx(elts_J0(2,:),:);
CJ0 = bndmesh.vtx(elts_J0(3,:),:);

ABCI0 = elts_I0;
ABCJ0 = elts_J0;

permI0 = repmat([1 2 3],size(elts_IJ0,1),1);
permJ0 = permI0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% COMMON VERTEX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[I1 J1] = find(Intmat == 1);

% Getting elements in matrix form corresponding to I1 and J1
elts_I1 = bndmesh.elt(I1,:);
elts_J1 = bndmesh.elt(J1,:);

% cell_elts_I1 = mat2cell(elts_I1,ones(size(I1,1),1));
% cell_elts_J1 = mat2cell(elts_J1,ones(size(J1,1),1));
% 
% % Apply intersection element wise
% intersection1 = cellfun(@intersect,cell_elts_I1,cell_elts_J1,'UniformOutput',false);
% % Apply setdiff element wise
% diffI1 = cellfun(@setdiff,cell_elts_I1,intersection1,'UniformOutput',false);
% diffJ1 = cellfun(@setdiff,cell_elts_J1,intersection1,'UniformOutput',false);
% 
% % Converting back to matrix
% intersection1 = cell2mat(intersection1);
% diffI1 = cell2mat(diffI1);
% diffI2 = cell2mat(diffI2);

[intersection1,diffI1,diffJ1] = rowWiseIntersectionDiff(elts_I1,elts_J1);

% Getting the vertices

% Common vertex! AI1 = AJ1
AI1 = bndmesh.vtx(intersection1,:);

BI1 = bndmesh.vtx(diffI1(:,1),:);
CI1 = bndmesh.vtx(diffI1(:,2),:);

BJ1 = bndmesh.vtx(diffJ1(:,1),:);
CJ1 = bndmesh.vtx(diffJ1(:,2),:);

ABCI1 = [intersection1 diffI1];
ABCJ1 = [intersection1 diffJ1];

% Finding permutation, perm is such that elt(perm) = ABC
permI1 = findPermVectorized(ABCI1,elts_I1);
permJ1 = findPermVectorized(ABCJ1,elts_J1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% COMMON EDGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[I2,J2] = find(Intmat == 2);

elts_I2 = bndmesh.elt(I2,:);
elts_J2 = bndmesh.elt(J2,:);

[intersection2,diffI2,diffJ2] = rowWiseIntersectionDiff(elts_I2,elts_J2);

% Getting the vertices

% Common edge! AI2 = AJ2, BI2 = BJ2
AI2 = bndmesh.vtx(intersection2(1,:),:);
BI2 = bndmesh.vtx(intersection2(2,:),:);

CI2 = bndmesh.vtx(diffI2,:);
CJ2 = bndmesh.vtx(diffJ2,:);

ABCI2 = [intersection2 diffI2];
ABCJ2 = [intersection2 diffJ2];

% Finding permutation, perm is such that elt(perm) = ABC
permI2 = findPermVectorized(ABCI2,elts_I2);
permJ2 = findPermVectorized(ABCJ2,elts_J2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% IDENTICAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[I3,J3] = find(Intmat == 3);

% Element I and J are identical
elts_IJ3 = bndmesh.elt(I3,:);

% Getting the vertices

% Identical element! AI3 = AJ3, BI3 = BJ3, CI3 = CJ3
AIJ3 = bndmesh.vtx(elts_IJ3(1,:),:);
BIJ3 = bndmesh.vtx(elts_IJ3(2,:),:);
CIJ3 = bndmesh.vtx(elts_IJ3(3,:),:);

ABCI3 = elts_IJ3;
ABCJ3 = elts_IJ3;

% Finding permutation, perm is such that elt(perm) = ABC
%permI3 = findPermVectorized(ABCI3,elts_IJ3);
%permJ3 = findPermVectorized(ABCJ3,elts_IJ3);
permI3 = repmat([1 2 3],size(elts_IJ3,1),1);
permJ3 = permI3;



