% Compute sphere sphere force at a distance using approximation for
% derivative
function out = sphsph_analytic_force(Qa,a,b,s)
    c = s+a+b;

    U = acosh((c^2-a^2-b^2)/(2*a*b));

    %caapr = Caa_pr(a,b,s);
    dCaads = Caa_pr_analytic(a,b,U,50);
    caa = Caa(a,b,U,50);

    out = -Qa^2*1/2/caa^2*dCaads/(4*pi);

end