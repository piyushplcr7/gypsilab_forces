function [M] = moments2(a,b,alpha,n)
%

ns = (0:(n-1))';
M = b.^(alpha + ns + 1)./(ns+alpha+1) - a.^(alpha + ns + 1)./(ns+alpha+1);



end

