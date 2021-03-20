function [ms] = nFoldMultiScreen(N,n)


Q = round(N/n);

originalPanel = mshSquare(Q,[1,2]);
singlePanel1 = originalPanel;
ang = 2*pi/n;
panels = cell(n,1);
for i = 1:n-1
    singlePanel2 = swap(rotate(singlePanel1,[-0.5,0,0],[0,1,0],ang));
    panels{i} = mshClean(union(singlePanel1,singlePanel2));
    singlePanel1 = swap(singlePanel2);
end
singlePanel2 = originalPanel;
if mod(n,2)==1
    singlePanel2 = swap(singlePanel2);
end
panels{end} = mshClean(union(singlePanel1,singlePanel2));


ms = multiScreen(panels);
ms = rotate(ms,[0,0,0],[1,0,0],pi/2);


end

