function[P] = mshEltVol2eltSurf(m,dm)

idx = 1:dm.nelt;
jdx = 0*idx;
val = 0*idx+1;

tol = 1e-12;
if size(m.elt,2)==3
    
    [~,I] = ismembertol(dm.vtx,m.vtx,'ByRows',true);
    b1 = ismember(m.elt(:,1),I);
    b2 = ismember(m.elt(:,2),I);
    b3 = ismember(m.elt(:,3),I);
    jdx = find(b1 + b2 + b3 >=2);
    
    
elseif size(m.elt,2)==4
    [~,I] = ismembertol(dm.vtx,m.vtx,'ByRows',true);
    b1 = ismember(m.elt(:,1),I);
    b2 = ismember(m.elt(:,2),I);
    b3 = ismember(m.elt(:,3),I);
    b4 = ismember(m.elt(:,4),I);
    jdx = find(b1 + b2 + b3 + b4 >=3);
else
    error('unavailable case.')
end

P = sparse(idx,jdx,val,dm.nelt,m.nelt);


end