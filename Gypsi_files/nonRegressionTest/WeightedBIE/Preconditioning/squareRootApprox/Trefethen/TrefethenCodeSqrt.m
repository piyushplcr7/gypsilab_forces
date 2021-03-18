A = pascal(5);
% change this for another matrix A
X = sqrtm(A);
% this cannot be changed
I = eye(size(A)); e = eig(A);
m = min(e); M = max(e);
% only for toy problems
k2 = m/M;
% elliptic functions parameter k^2
Kp = ellipke(1-k2);
for N = 5:5:20
    t = 1i*(.5:N)*Kp/N;
    [sn, cn, dn] = ellipj(imag(t),1-k2);
    cn = 1./cn; dn = dn.*cn; sn = 1i*sn.*cn;
    w = sqrt(m)*sn;
    dzdt = cn.*dn;
    S = zeros(size(A));
    for j = 1:N
        S = S - inv(A-w(j)^2*I)*dzdt(j);
    end
    S = (-2*Kp*sqrt(m)/(pi*N))*A*S;
    error(N) = norm(S-X)/norm(X);
    fprintf('%4d %10.2e\n', N, error(N))
end