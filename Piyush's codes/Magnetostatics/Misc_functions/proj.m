% To project a function a to space, its value at the quadrature points must
% be passed in A.
function coeffs = proj(A,Gamma,space)
    % Mass matrix
    M = integral(Gamma,space,space);

    uqmat = space.uqm(Gamma);
    [X,W] = Gamma.qud;

    assert(size(X,1)==size(A,1));

    if (iscell(uqmat))
        assert(size(A,2)==3);
        % Computing the RHS, ie (A,b_i)
        rhs = sum(W.*(A(:,1).*uqmat{1} + A(:,2).*uqmat{2} + A(:,3).*uqmat{3}),1)';
    else
        assert(size(A,2)==1)
        rhs = sum(W.*(A.*uqmat),1)';
    end

    coeffs = M\rhs;
end