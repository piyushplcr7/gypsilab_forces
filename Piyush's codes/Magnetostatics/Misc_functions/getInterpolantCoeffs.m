% Interpolation function for RWG or NED spaces

function coeffs = getInterpolantCoeffs(f,space)
    type = space.typ;
    dofs = space.dof; % edge centers
    bndmesh = space.msh;
    [edgmesh,elt2edg] = mshEdge(bndmesh);
    tgts = edgmesh.tgt;

    nrms = 0 * tgts;

    

    edgnrms = bndmesh.nrmEdg;
    edgtgts = bndmesh.tgt;
    coeffs = 
    switch type
        case "RWG"
            plotTgt(bndmesh)

        case "NED"

    end



end