function [ phi ] = DarbasPadeSqrt( u,Np,theta,I,D)

% Approximation of R = I * sqrtm(I^{-1}D)

[ C0,Aj,Bj ] = rotatingPadeRacine(Np,theta);
phi = C0*I*u;

for j = 1:Np
    phi_j = (Bj(j)*(D-I) + I)\(I*u);
    phi = phi + Aj(j)*(D-I)*phi_j;
end


end
