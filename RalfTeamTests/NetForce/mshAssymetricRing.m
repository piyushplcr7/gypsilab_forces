function mesh = mshAssymetricRing(N,rad1,rad2)


rad = rad2-rad1;
% Radial discretisation
dr = sqrt(pi*rad^2/N);
dr = rad/ceil(rad/dr);
r  = rad1:dr:rad2;

% Angular uniform discretization
rho = cell(length(r),1); theta = rho;
for ir = 1:length(r)
    alpha = (r(ir) -rad1/2)/(rad2-rad1/2);
    dtheta = dr/r(ir);
    dtheta = 2*pi/ceil(2*pi/dtheta);    
    theta{ir} = (0:dtheta:2*pi-dtheta)';
    rho{ir}   = r(ir).*(1 + alpha*(1 + 0.5*cos(theta{ir})));
end

% Carthesian coordinates
[x,y] = pol2cart(cell2mat(theta),cell2mat(rho));
X     = [0,0; x y];

% Unicity test
tmp = unique(X,'rows','stable');
if (max(abs(X-tmp)) > 1e-12)
    error('mshDisk : non unicity of vertices')
end
   
% Delaunay triangulation
DT = delaunayTriangulation(X(:,1),X(:,2));

% Final mesh
elt  = DT.ConnectivityList;
vtx  = [DT.Points,zeros(size(DT.Points,1),1)];
mesh = msh(vtx,elt);

I = find(mesh.vtx(:,1).^2 + mesh.vtx(:,2).^2 < rad1.^2/2);
ind = and(and(~ismember(mesh.elt(:,1),I),~ismember(mesh.elt(:,2),I)),...
    ~ismember(mesh.elt(:,3),I));
mesh = mesh.sub(ind);

end
