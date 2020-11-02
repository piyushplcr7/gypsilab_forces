function [ v ] = TrefethenSqrtGalerkin(Ah,N,u,Ih,min_m,max_M)
% returns v an approximation of sqrtm(A)
% or v approx sqrtm(A)*u
% if u is not empty. 

k2 = min_m/max_M; % elliptic functions parameter k^2
Kp = ellipke(1-k2);
t = 1i*(.5:N)*Kp/N;
[sn, cn, dn] = ellipj(imag(t),1-k2);
cn = 1./cn; dn = dn.*cn; sn = 1i*sn.*cn;
w = sqrt(min_m)*sn;
dzdt = cn.*dn;

if ~isempty(u)
    
    Sv = zeros(size(Ah,1),1);
    Ahu = Ah*u;
    for j = 1:N
        Sv = Sv - Ih*((Ah-w(j)^2*Ih)\Ahu)*dzdt(j);
    end
    v = (-2*Kp*sqrt(min_m)/(pi*N))*Sv;
else
    
    S = zeros(size(Ah));
    for j = 1:N
        S = S - Ih*inv(Ah-w(j)^2*Ih)*dzdt(j);
    end
    v = (-2*Kp*sqrt(min_m)/(pi*N))*S*Ah;
    
end

