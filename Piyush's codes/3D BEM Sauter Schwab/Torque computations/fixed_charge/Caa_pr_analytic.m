% Evaluate Caa_pr analytically for identical spheres
function out = Caa_pr_analytic(a,b,U,n)
    out1 = 0;
    out2 = 0;
    %c = sqrt(2*(cosh(U)+1))*a;
    c = sqrt(2*a*b*cosh(U)+a^2+b^2);

    dUds = c/a/b/sinh(U);
    for i = 0:n
        out1 = out1 + 1/( a*sinh(i*U) + b*sinh((i+1)*U) );
        out2 = out2 + 1/( a*sinh(i*U) + b*sinh((i+1)*U) )^2 * (a * cosh(i*U) *i + b* cosh((i+1)*U) * (i+1));

    end

    out1 = out1 * a * b * cosh(U) * dUds;

    out2 = -out2 * a *b * sinh(U) * dUds;

    out = out1 + out2;

end