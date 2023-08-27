% integral of sin(a1x) sin(a2x) from -Npi to Npi
function out = Ja1a2(a1,a2,N)
    if a1 ~= a2
        out = 0;
        return;
    elseif a1 == a2 && a1>0
        out = N * pi;
        return;
    else % a1 == a2 == 0
        out = 0;
    end
end