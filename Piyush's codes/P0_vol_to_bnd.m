function tr_opr = P0_vol_to_bnd(mesh)
%% Creating projection matrix linking boundary triangles to boundary edges
bnd_mesh = mesh.bnd;
% Number of elements in boundary and volume mesh
Nelt_b = size(bnd_mesh.elt,1);
Nelt_v = size(mesh.elt,1);
tr_opr = sparse(Nelt_b, Nelt_v);

% Looping over all boundary elements and finding the parent volume element
for i = 1 : Nelt_b
    bvindx = bnd_mesh.elt(i,:);
    % Nodes of the boundary element
    bvpts = bnd_mesh.vtx(bvindx,:);
    % Checking which volume mesh nodes match with selected boundary nodes
    distpt1 = vecnorm(mesh.vtx-ones(size(mesh.vtx,1),1).*bvpts(1,:),2,2);
    distpt2 = vecnorm(mesh.vtx-ones(size(mesh.vtx,1),1).*bvpts(2,:),2,2);
    % Indices for matching nodes
    j1 = find(distpt1 == 0);
    j2 = find(distpt2 == 0);
    % Checking which volume element contains both nodes
    for j = 1:Nelt_v
       cur_elt = mesh.elt(j,:);
       % Checking if coincides
       check1 = cur_elt(1) == j1 | cur_elt(2) == j1 | cur_elt(3) == j1;
       check2 = cur_elt(1) == j2 | cur_elt(2) ==j2 | cur_elt(3) == j2;
       if check1 & check2
           disp('Found match');
          tr_opr(i,j)=1; 
       end
    end
end

end