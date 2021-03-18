function[M] = moments(d,n,a,b)

% Computes the integrals \int_{a}^b u/(u^2 + d^2) u^i du
% for i = 1..n
% using the fact that 2u/(u^2 + d^2) = 1/(u + id) + 1/(u-id)
% and the recurrence relation u^n/(u + id) = u^(n-1) -id u^(n-1)/(u+id)


M1 = zeros(length(d),n);
M2 = zeros(length(d),n);
BN = zeros(length(d),n);
AN = zeros(length(d),n);
M1(:,1)= log((b+1i*d)./(a+1i*d));
M2(:,1)= log((b-1i*d)./(a-1i*d));
AN(:,1)= 1;
BN(:,1)= 1;

for i = 2:n
    AN(:,i) = AN(:,i-1)*a;
    BN(:,i) = BN(:,i-1)*b;
    M1(:,i) = BN(:,i)/(i-1) - AN(:,i)/(i-1) - 1i*d.*M1(:,i-1);
    M2(:,i) = BN(:,i)/(i-1) - AN(:,i)/(i-1) + 1i*d.*M2(:,i-1);
end


M = (M1 + M2)/2;


end