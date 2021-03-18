function [] = pltNaturalNrm(varargin)

m = varargin{1};
spc  = 'r';
if (nargin == 2)
    spc = varargin{2};
end


fce = m.fce;
fce.elt = sort(fce.elt,2);
Xctr = fce.ctr;
Vnrm = cross(fce.vtx(fce.elt(:,2),:) - fce.vtx(fce.elt(:,1),:),...
    fce.vtx(fce.elt(:,3),:) - fce.vtx(fce.elt(:,1),:));
Vnrm = Vnrm./sqrt((Vnrm(:,1).^2 + Vnrm(:,2).^2+Vnrm(:,3).^2 ));


quiver3(Xctr(:,1),Xctr(:,2),Xctr(:,3),Vnrm(:,1),Vnrm(:,2),Vnrm(:,3),...
    'Color',spc);




end

