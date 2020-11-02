function [logR,rlogR,gradlogR,sX,d,omegaX] = WdomSemiAnalyticInt2D(Ydom,X,A,B,eltNum)
% Exact integration of I = \int_[a,b] Gd((s(x)-s(y))/s'(x)) s'(y) dy
% where 1/omega(X) = s'(x) is the weight function on Ydom and
% X = A + x t + d n where t and n are the tangent and normal vectors on
% [A,B].
% Gd is 
% - 1/2 ln(d^2 + x^2)       (logR)
% - d/2 ln(d^2 + x^2)       (rlogR, normal part)
% - x/2 ln(d^2 + x^2)   (rlogR, tangential part)
% - d/(d^2 + x^2)           (gradlogR, normal part)
% - x/(d^2 + x^2)         (gradlogR, tangential part)
% 
% Knowing the primitive F of Gd, the integral is explicitly given by 
% I = s'(x) * {F[s'(x)^(-1)(s(b) - s(x))] - F[s'(x)^(-1)(s(a) - s(x))]}



Tel = (B-A)./norm3D(B-A);
Nel = [Tel(2) -Tel(1) 0];


sA = Ydom.sA(eltNum);
sB = Ydom.sB(eltNum);
sX = Ydom.r_1(eltNum,X);
omega = Ydom.w;
omegaX = omega(X);


% Integral parametrization
d = abs(scal3D(X - A,Nel));


% logR
logR = omegaX.*(F1((sB - sX)./omegaX) - F1((sA - sX)./omegaX));

% rlogR
rlogRn = d.*logR;
rlogRt = omegaX.*(F2((sB - sX)./omegaX) - F2((sA - sX)./omegaX));

rlogR(:,1) = rlogRn*Nel(1) + rlogRt*Tel(1);
rlogR(:,2) = rlogRn*Nel(2) + rlogRt*Tel(2);
rlogR(:,3) = rlogRn*Nel(3) + rlogRt*Tel(3);

% gradlogR (R/r^2)

gradlogRn = omegaX.*(F3((sB - sX)./omegaX) - F3((sA - sX)./omegaX));
gradlogRt = omegaX.*(F4((sB - sX)./omegaX) - F4((sA - sX)./omegaX));

gradlogR = zeros(size(X));
gradlogR(:,1) = gradlogRn*Nel(1) + gradlogRt*Tel(1);
gradlogR(:,2) = gradlogRn*Nel(2) + gradlogRt*Tel(2);
gradlogR(:,3) = gradlogRn*Nel(3) + gradlogRt*Tel(3);





function y = F1(x)

y = 1/2*x.*log(x.^2+d.^2) - x + d.*atan(x./d);
y(d==0) = xlog(x(d==0)) - x(d==0);

end


function y  = F2(x)
% Primitive of x/2 ln(x^2 + d^2)

y = 1/4*(xlog(d.^2 + x.^2) - x.^2);

end

function y = F3(x)
% Primitive of d/(d^2 + x^2)
y = atan(x./d);
y(abs(d) < 1e-13) = 0;

end

function y = F4(x)
% Primitive of x/(d^2 + x^2)
y = 1/2*log(x.^2 + d.^2);

end

function y = xlog(x)
y    = x.*log(abs(x));
y(abs(x)<1e-13) = 0;
end


function[res] =  I1(alpha,b,d,x)
% Returns \int_{0}^b y^(-alpha) ln(\sqrt(d^2 + (x-y)^2) dy

TB = 1/(1-alpha)*b^(1-alpha)/2*(log((b-x).^2 + d.^2));
z1 = x + 1i*d;
z2 = x - 1i*d;

res = TB - 1/(1-alpha)*(I4(1-alpha,b,z1) + I4(1-alpha,b,z2))/2;

end

function[res] = I2(alpha,b,d,x)
% Returns \int_{0}^b  y^(-alpha)d / ((y - x)^2 + d^2) dy
% d/((y-x)^2 + d^2) = 1/(2i)*(1/(y-x+id) - 1/(y-x - id))

res = 0*d;
res(d==0) = 0;
dnz = d(d~=0);
xnz = x(d~=0);
z1 = xnz - 1i*dnz;
z2 = xnz + 1i*dnz;
res(d~=0) = 1/(2i)*(I4(-alpha,b,z1) - I4(-alpha,b,z2));

end

function[res] = I3(alpha,b,d,x)
% Returns \int_{0}^b y^(-alpha)(y-x)/ ((y - x)^2 + d^2) dy
% (y-x)/((y-x)^2 + d^2) = 1/2*(1/(y-x-id) + 1/(y-x + id))

z1 = x + 1i*d;
z2 = x - 1i*d;
res = 1/2*(I4(-alpha,b,z1) + I4(-alpha,b,z2));
end

function[res] = I4(alpha,b,z)
% Returns \int_{0}^b y^(alpha)/(y - z) dy
Naccu = 5;
res = -(b^(alpha+1))./((alpha+1)*z).*hypergeom([1.,alpha+1.],alpha + 2.,b./z);
res(imag(z) == 0) = specialIntegral1(-alpha,b,z(imag(z)==0),Naccu);

end






% function[res] = aux1(y,d)
% % Primitive of y-> ln(sqrt(y^2 + d^2))
% res = -y + y.*log(sqrt(y.^2 + d.^2) + d.*atan(y./d));
% res(abs(d)<1e-12) = -y + ylog(y);
% end
% function[res] = aux2(y)
% % Primitive of y-> y * ln(|y|)
% res = 1/4*(y.*ylog(y.^2) - y.^2);
% 
% end
% 
% function[res] = aux3(t,u,Y,alpha,beta,b,d,W,Ydom)
% % regularizes 1/2*\int_0^b y^-beta \ln |d^2 + (x-y)^2|
% % alpha is the power singularity of the Wdom object.
% 
% k = @(u)(1/(1 - beta)*abs(u).^(1-beta));
% dk = @(u)(abs(u).^(-beta));
% H = @(u,d)(dk(t).*(aux1((k(t) - k(u))./dk(t),d)));
% I0exact = -H(b,d) + H(0,d);
% 
% T = t*ones(size(u'));
% U = ones(size(t))*u';
% D = d*ones(size(u'));
% TU = (k(T) - k(U))./dk(T);
% fTU = 1/2*log(D.^2 + TU.^2);
% fTU(or(isinf(fTU),isnan(fTU))) = 0;
% I0approx = fTU*(W./Ydom.w(Y).*(abs(u).^(-beta)));
% 
% 
% res = I0exact - I0approx;
% 
% end
% 
% function[res] = aux4(t,u,Y,alpha,b,d,W,Ydom)
% % regularizes 1/2*\int_0^b y^-alpha * 1/(y-x)
% 
% 
% Prim = @(u)(-t.^(-alpha).*betainc(u./t,1 - alpha,0));
% I0exact = Prim(b) - Prim(0.0);
% T = t*ones(size(u'));
% U = ones(size(t))*u';
% D = d*ones(size(u'));
% TU = (k(T) - k(U))./dk(T);
% fTU = -1./TU;
% fTU(or(isinf(fTU),isnan(fTU))) = 0;
% I0approx = fTU*(W./Ydom.w(Y).*(abs(u).^(-alpha)));
% 
% 
% res = I0exact - I0approx;
% 
% end
% 
% 
% function[res] = ylog(y)
% res = y.*log(abs(y));
% res(isnan(res)) = 0;
% res(isinf(res)) = 0;
% end



end













