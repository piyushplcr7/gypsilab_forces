function [gammaU,Wh] = trace_P0P1(Vh)
% Computes the restriction of uh to the mesh m.

m = Vh.msh;
dm = m.bnd;

if strcmp(Vh.typ,'P0')
    P = mshEltVol2eltSurf(m,dm);
    Wh = fem(dm,'P0');
    if nargin==2
        gammaU = P*uh;
    else
        gammaU = P;
    end
elseif strcmp(Vh.typ,'P1')
    Wh = fem(dm,'P1');
    [~,I1,I2] = intersect(Wh.unk,Vh.unk,'rows','stable');
    P = sparse(I1,I2,1,length(Wh),length(Vh));
    gammaU = P;    
else
    error('Unavailable case');
end


end

