function [ phi ] = padePrecondDarbas( l,Np,theta,keps,I,D)


k = real(keps);
[ C0,Aj,Bj ] = rotatingPadeRacine(Np,theta);

if isempty(l)
    phi = 1i*k*C0*I;
    for j = 1:Np
        Phi_j = (-Bj(j)/keps^2*D + I)\I;
        phi = phi + (-1i*k*Aj(j)/keps^2)*(D*Phi_j);
    end
else
    phi = 1i*k*C0*I*l;
    
    for j = 1:Np
        phi_j = (-Bj(j)/keps^2*D + I)\(I*l);
        phi = phi + (-1i*k*Aj(j)/keps^2)*(D*phi_j);
    end
end



end

