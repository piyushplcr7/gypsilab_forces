% Function to compute Caa'
function out = Caa_pr(a,b,s)
    c = s+a+b;
    delta = 1e-12;
    cp = s+a+b+delta;
    U = (c^2-a^2-b^2)/(2*a*b);
    Up = (cp^2-a^2-b^2)/(2*a*b);

    val = Caa(a,b,U,50);
    valp = Caa(a,b,Up,50);
    out = (valp-val)/delta;

end