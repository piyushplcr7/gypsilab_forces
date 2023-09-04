% Computation done in this function is futile becaus ethe shape derivatives
% are not expected to match for velocity fields that don't create an
% isometry. They are expected to be true only when the fileds are constant
% fields or rotational fields.

function [] = compareSPVPBEM(fnameVP,fnameSP,trim)
    load(fnameVP);
    hvalsVP = hvals(trim);
    shape_derivatives_bem_VP = shape_derivatives_bem(trim,:);

    load(fnameSP);
    hvalsSP = hvals(trim);
    shape_derivatives_bem_SP = shape_derivatives_bem(trim,:);

    % dual norm of the error between BEM SP and BEM VP
    gramat = gramat_abc_alpha();
    diagonal = diag(gramat)';
    diagonal = 1./diagonal;

    errmat = shape_derivatives_bem_VP - shape_derivatives_bem_SP;
    
    Ginv_errs = diagonal.*errmat;

    dualnormerr = sqrt(dot(errmat,Ginv_errs,2));

    figure;
    loglog(hvals,dualnormerr);


end