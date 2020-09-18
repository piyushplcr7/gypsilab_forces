function[m] = meshCurve(c,N,varargin)
% Create a mesh of a curve c with N elts
%
% The mesh density can be controlled in two ways:
% - The user can provide a change of variables z together with its interval
% of definition I_z. The mesh will be generated as a uniform mesh for the
% curve c parametrized on I_z by u-> (c.x(z(u)), c.y(z(u)), if z([I_z]) =
% c.I. In case z([I_z]) != c.I, the curve c is first affinely
% reparametrized to ensure this condition.
% - Alternatively, the user can provide a density function d together with
% an interval of definition I_d. A change of variables z defined on
% I_z = [0, q(I_d(2))] is then computed from this density by the formula
% z = q^{-1} where q(t) = int_{I_d(1)}^t d(u) du.
% To perform this computation, the density function is approximated
% by its linear interpolation at sample points which by default are chosen
% as linspace(I_d(1),I_d(2),M) where M = 2*N, and q is in turn approximated
% by its linear interpolant at those same grid points.

p = inputParser;
p.addOptional('density',{[],[]});
p.addOptional('varChange',{[],[]});
p.parse(varargin{:});




if ~isempty(p.Results.density{1})
    d = p.Results.density{1}; I_d = p.Results.density{2};
    a = I_d(1); b = I_d(2);
    samples = linspace(a,b,2*N);
    [F_t,t] = repartitionFunction(samples,d(samples));
    z = @(s)(interp1(F_t,t,s));
    I_z = [F_t(1),F_t(end)];
elseif ~isempty(p.Results.varChange{1})
    z = p.Results.varChange{1};
    I_z = p.Results.varChange{2};
end


if exist('z','var')&&~isempty(z)
    alpha = I_z(1);
    beta = I_z(2);
    a = min(z(alpha),z(beta));
    b = max(z(alpha),z(beta));
    bool1 = abs(a-c.I(1))<1e-10;
    bool2 = abs(b-c.I(2))<1e-10;
    if ~and(bool1,bool2)
        c = c.affineReparam([a,b]);
    end
    if c.closed
        t = z(linspace(alpha,beta,N));
        elt = [(1:N)', [(2:N)';1]];
    else
        t = z(linspace(alpha,beta,N+1));
        elt = [(1:N)',(2:N+1)'];
    end
else
    if c.closed
        t = linspace(c.I(1),c.I(2),N+1);
        t = t(1:end-1);
        elt = [(1:N)', [(2:N)';1]];
    else
        t = linspace(c.I(1),c.I(2),N+1);
        elt = [(1:N)',(2:N+1)'];
    end
end

t = t(:);
x = c.x(t);
y = c.y(t);
vtx = [x(:), y(:)];
m = msh(vtx,elt);

end

function [F_t,t] = repartitionFunction(t,d_t)
% Computes the values F_t(i) = F(t(i)) of the function F(s) defined by 
% F(s) = \int_{t(1)}^s d(u) du 
% where d(u) is the piecewise linear function defined such that f(t(i)) =
% d_t(i). 

input_shape = size(t);
t = t(:);
d_t = d_t(:);


A = t(1:end-1);
B = t(2:end);
Alpha = d_t(1:end-1);
Beta = d_t(2:end);

val = (Alpha + Beta).*(B - A)/2;

F_t = [0; cumsum(val)];

F_t = reshape(F_t,input_shape);
t = reshape(t,input_shape);

end



