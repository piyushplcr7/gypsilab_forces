function [wA,wB,wC,ang] = standardTri(zA,zB,zC)
% Returns the complex coordinates of the verties of a triangle 
% not containing the origin obtained form (zA,zB,zC) by a rotation, and 
% reordering the points such that:
% the (-pi,pi] argument of wA, wB and wC, phiA,phiB and phiC, are such that
% 0 = phiA <= phiB <=phiC

phiA = arg(zA); phiB = arg(zB); phiC = arg(zC);
rTri = [abs(zA),abs(zB),abs(zC)];

phi = [phiA,phiB,phiC];
ang = -phi(1);
[phi,J] = sort(myMod2pi(phi - phi(1))); % rotation of angle -phi(1)
rTri = rTri(J);

ang = ang - phi(1);
phi = myMod2pi(phi - phi(1));

wA = rTri(1);
wB = rTri(2)*exp(1i*phi(2));
wC = rTri(3)*exp(1i*phi(3));

end

