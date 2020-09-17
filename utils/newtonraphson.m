function [x] = newtonraphson(f,x0,df,varargin)

% Looks for a zero of f near x0 using the Newton-Raphson method algorithm

MAXIT = 20;
TOL = 1e-6;

p = inputParser;
p.addOptional('maxit',MAXIT);
p.addOptional('tol',TOL);

p.parse(varargin{:});


x = x0;
for i = 1:p.Results.maxit
    crit = abs(f(x)./df(x));
    x = x - f(x)./df(x);
    if crit<p.Results.tol
        return
    end
end

warning(['Newton-Raphson algorithm did' ...
        ' not converge after %s iterations. Crit = %s']...
    ,i,crit);

end

