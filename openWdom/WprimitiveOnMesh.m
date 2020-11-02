function [p] = WprimitiveOnMesh(m,y,omega,singVtx,singPow,N)
% For each y(i) on the semgent [A(i),B(i)], computes p(y) = r^{-1}(y(i)) -
% r^{-1}(A(i))
% where r^{-1}(y) = \int_{0}^{s^{-1}(y)} omega(s) ds
% and where s is an arclength on m (if m is an open curve, then s(0) must
% be one of the two edges.
% Thus, since [A(i),B(i)] is a segment, p(y) = \int_{[A,y]} omega(x) dx

% Sing contains locations and power of singularities in omega.


computeMatlabAnswer = false;


M = size(y,1);
[A,B] = m.ABCD;
n = m.nrm;
t = m.tgt;

Q = M/m.nelt;
As = repeatLines(A,Q);
Bs = repeatLines(B,Q);
ns = repeatLines(n,Q);
ts = repeatLines(t,Q);
%
% % Matlab answers
if computeMatlabAnswer
    for i = 1:M
        f = @(s)(reshape(omega(As(i,:) + s(:)*ts(i,:)),size(s,1),size(s,2)));
        yi = scal3D(y(i,:)-As(i,:),ts(i,:));
        p2(i,1) = integral(f,0,yi);
    end
    
end


h = m.stp;
[I,r] = rangesearch(y,singVtx,3*h(2));
% Points that are far from singularities

[x,w] = domGauss_Legendre1D(N,0,1);
x = sort(x);

X = repeatLines(As,N);

for i = 1:3
    X(:,i) = X(:,i) + reshape(((y(:,i) - As(:,i))*x')',1,N*M)';
end

w = kron(norm3D(y - As),w');
omegaX = reshape(omega(X),N,M)';

p = sum(omegaX.*w,2);

for j = 1:size(singVtx(:,1))
    Sj = singVtx(j,:);
    alphaj = singPow(j);
    if ~isempty(I{j})
        for ind = 1:length(I{j})
            i = I{j}(ind);
            Ai = As(i,:);
            ni = ns(i,:);
            ti = ts(i,:);
            di = scal3D(ni,Sj - Ai);
            si = scal3D(ti,Sj - Ai);
            Li = scal3D(ti,y(i,:)-Ai); % always positive
            if r{j}(ind)<8*Li
                
                if abs(di) > 1e-12
                    if and(si + di>=0,si-di <= Li)
                        thetaA = max(-pi/4,atan(-si/di));
                        thetaB = min(pi/4,atan((Li-si)/di));
                        Ni = fix(1.5*N*abs(alphaj*log(di)))+1;
                        [the,w] = Gauss_Legendre1D(Ni,thetaA,thetaB);
                        the = sort(the);
                        X = zeros(Ni,3);
                        for k = 1:3
                            X(:,k) = Sj(k) - di*ni(k) + di*tan(the)*ti(k);
                        end
                        W = di*(1 + tan(the).^2).*w;
                        I1 = sum(omega(X).*W);
                        if atan((Li-si)/di) > pi/4
                            Ni = abs(fix(log(di/(Li - si))*N)+1);
                            [x,w] = Gauss_Legendre1D(Ni,log(di),log(Li-si));
                            x = sort(x);
                            X = zeros(Ni,3);
                            for k = 1:3
                                X(:,k) = Sj(k) - di*ni(k) + exp(x)*ti(k);
                            end
                            w = exp(x).*w;
                            I2 = sum(w.*omega(X));
                        else
                            I2 = 0;
                        end
                        if atan(-si/di) < -pi/4
                            [x,w] = Gauss_Legendre1D(20,log(di),log(si));
                            X = zeros(20,3);
                            for k = 1:3
                                X(:,k) = Sj(k) - di*ni(k) - (exp(x))*ti(k);
                            end
                            w = exp(x).*w;
                            I3 = sum(w.*omega(X));
                        else
                            I3 = 0;
                        end
                        p(i) = I1 + I2 + I3;
                    end
                else
                    % Sj is positioned at 0.
                    
                    a = -si;
                    b = Li - si;
                    N1 = max(N,fix(N*abs(a)/Li)+1);
                    N2 = max(N,fix(N*abs(b)/Li)+1);
                    [x1,w1] = gaussQuadPower(-alphaj,N1,a);
                    [x2,w2] = gaussQuadPower(-alphaj,N2,b);
                    X1 = zeros(N1,3);
                    X2 = zeros(N2,3);
                    for k = 1:3
                        X1(:,k) = Sj(k) + x1*ti(k);
                        X2(:,k) = Sj(k) + x2*ti(k);
                    end
                    sa = (a<0)*2 - 1;
                    sb = (b>0)*2 - 1;
                    omegaX1 = omega(X1).*(abs(x1).^(-alphaj));
                    omegaX1(isnan(omegaX1))=0;
                    omegaX2 = omega(X2).*(abs(x2).^(-alphaj));
                    omegaX2(isnan(omegaX2))=0;
                    p(i) = sa*sum(omegaX1.*w1) ...
                        + sb*sum(omegaX2.*w2);
                    
                    if computeMatlabAnswer
                        f = @(s)(reshape(omega(As(i,:) + s(:)*ts(i,:)),size(s,1),size(s,2)));
                        if abs(si) < 1e-13
                            p2(i,1) = integral(f,si,Li);
                        elseif abs(si - Li) < 1e-13
                            p2(i,1) = integral(f,0,si);
                        else
                            p2(i,1) = integral(f,0,si) + integral(f,si,Li);
                        end
                    end
                end
            end
            
        end
        
    end
end
% %
if computeMatlabAnswer
    disp(max(abs(p-p2)));
end
end

