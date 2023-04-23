% Computing the capacitance Caa
function out = Caa(a,b,U,n)
    out = 0;
    for i = 0:n
        out = out + 1/(a*sinh(i*U)+b*sinh((i+1)*U));
    end
    out = out * a*b*sinh(U);
end