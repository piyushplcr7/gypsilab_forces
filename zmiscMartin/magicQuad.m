function [z,w] = magicQuad(zA,zB,zC,NgaussPhi,NgaussR)

z = zeros(2*NgaussPhi,NgaussR);
w = zeros(2*NgaussPhi,NgaussR);

[wA,wB,wC,ang] = standardTri(zA,zB,zC);
% plot([wA wB wC wA],'k-');
phiB = arg(wB);

phiC = arg(wC);
[phi1,wphi1] = domGauss_Legendre1D(NgaussPhi,0,phiB);
[phi2,wphi2] = domGauss_Legendre1D(NgaussPhi,phiB,phiC);
phi = [phi1;phi2];
wphi = [wphi1;wphi2];

xA = real(wA); yA = imag(wA);
xB = real(wB); yB = imag(wB);
xC = real(wC); yC = imag(wC);
for j = 1:(2*NgaussPhi)
    if phi(j)< phiB % First sector
        tb =@(u)( -(yA*cot(u) - xA)./((yB-yA)*cot(u) - (xB - xA)));
        wAB = wA + tb(phi(j))*(wB - wA);
%         plot(wAB,'r*');
        rAB = abs(wAB);
        
        tc = @(u)(-(yA*cot(u) - xA)./((yC-yA)*cot(u) - (xC - xA)));
        wAC = wA + tc(phi(j))*(wC - wA);
%         plot(wAC,'r*');
        rAC = abs(wAC);
        r1 = min(rAB,rAC); r2 = max(rAB,rAC);
    else
        tb = -(yC*cot(phi(j)) - xC)/((yB-yC)*cot(phi(j)) - (xB - xC));
        wCB = wC + tb*(wB - wC);
%         plot(wCB,'r*');
        rCB = abs(wCB);
        
        ta = -(yC*cot(phi(j)) - xC)/((yA-yC)*cot(phi(j)) - (xA - xC));
        wCA = wC + ta*(wA - wC);
%         plot(wCA,'r*');
        rCA = abs(wCA);
        r1 = min(rCB,rCA); r2 = max(rCB,rCA);
    end
    [thet,wr] = domGauss_Legendre1D(NgaussR,acos(r2),acos(r1));
    for k = 1:NgaussR
        z(j,k) = cos(thet(k))*exp(1i*phi(j));
        w(j,k) = wphi(j)*wr(k)*cos(thet(k));
%         plot(z(j,k),'g*');
    end
    
    
    
end


z = z*exp(-1i*ang);
% plot(z,'b*');







end

