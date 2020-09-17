function [r] = compositeIntegral(f,a,b,varargin)
% Compute the integral of f using a composite Gaussian quadrature rule on [a,b]

p = inputParser;
p.addOptional('N',100);
p.addOptional('gss',3);
p.addOptional('tol',1e-6);
p.parse(varargin{:});
gss = p.Results.gss; N = p.Results.N;

t = linspace(a,b,N+1);
vtx = [t(:), 0*t(:), 0*t(:)];
elt = [1:N;2:(N+1)]';
m = msh(vtx,elt);
d = dom(m,gss);

F = @(Z)(f(Z(:,1)));

r = integral(d,F);


end

