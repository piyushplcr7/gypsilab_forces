function [x] = invFunc(F_dF,y,x0,tol)

if ~exist('tol','var')||isempty(tol)
    tol = 1e-3;
end

x = x0;
[fx,dfx] = F_dF(x);

err = max(abs(fx-y));

while err> tol
    x = x - (fx - y)./dfx;
    [fx,dfx] = F_dF(x);
    err = max(abs(fx-y));
end




end

