function [Q] = multPoly(A,B)


n = length(A)-1; % P1 of degree n
m = length(B)-1; % P2 of degree m

Q = zeros(1,n+m+1); % Q of degree n+m
% Constant term on the left.

for k = 0:(n+m)
    ind = k - (0:(n)) + 1;
    nz = and(ind>0,ind<=m+1);
    
    btilde = zeros(1,n+1);
    
    btilde(nz) = B(ind(nz));
    Q(k + 1) = sum(A.*btilde);
    
    
    
end




end

