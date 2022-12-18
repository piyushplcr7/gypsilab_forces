function y = reconstruct(coeffs,Gamma,space)
    uqmat = space.uqm(Gamma);

    if (iscell(uqmat))
        y1 = sum(coeffs'.*uqmat{1},2);
        y2 = sum(coeffs'.*uqmat{2},2);
        y3 = sum(coeffs'.*uqmat{3},2);
        y = [y1 y2 y3];
    else
        y = sum(coeffs'.*uqmat,2);
    end
end