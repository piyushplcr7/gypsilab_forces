function [nxU,Wh] = nx_RWGNED(Vh)

assert(size(Vh.msh.elt,2)==3);

if strcmp(Vh.typ,'NED')
    Wh = fem(Vh.msh,'RWG');
    nxU = fem_eye(Vh);
else
    Wh = fem(Vh.msh,'NED');
    nxU = -fem_eye(Vh);
end


end

