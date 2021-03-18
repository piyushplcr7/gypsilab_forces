function [gradU,Wh] = grad_P1(Vh)

m = Vh.msh;
Wh = fem(m,'NED');
[ed,~] = m.edg;
idx = (1:Wh.ndof)';
gradU = sparse(idx,ed.elt(:,1),1,Wh.ndof,Vh.ndof)...
    -sparse(idx,ed.elt(:,2),1,Wh.ndof,Vh.ndof) ;

if size(m.elt,2)==4
    gradU = -gradU;
end

end

