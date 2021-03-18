function [ v ] = TrefethenSqrt(A,N,u,min_m,max_M)
% returns v an approximation of sqrtm(A)
% or v approx sqrtm(A)*u
% if u is not empty. 

I = eye(size(A));
k2 = min_m/max_M; % elliptic functions parameter k^2
Kp = ellipke(1-k2);
t = 1i*(.5:N)*Kp/N;
[sn, cn, dn] = ellipj(imag(t),1-k2);
cn = 1./cn; dn = dn.*cn; sn = 1i*sn.*cn;
w = sqrt(min_m)*sn;
dzdt = cn.*dn;

if ~isempty(u)
    
    v = zeros(size(A,1),1);
    for j = 1:N
        v = v - dzdt(j)*((A-w(j)^2*I)\u);
    end
    v = (-2*Kp*sqrt(min_m)/(pi*N))*(A*v);
else
    
    S = zeros(size(A));
    for j = 1:N
        S = S - inv(A-w(j)^2*I)*dzdt(j);
    end
    v = (-2*Kp*sqrt(min_m)/(pi*N))*A*S;
    
end

