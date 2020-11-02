function [I] = specialIntegral1(alpha,b,X,N)


% Returns an accurate approximation of the integral
% int_{0}^b y^{-alpha}/(y-x) dy
% where alpha < 1
% x > 0, x ~= b, 
% understood in principal value if x lies in (0,b)
% N is a parameter controlling accuracy. 


% For now, consider that X = X(:)
X = X(:);

I1 = zeros(size(X));
I2 = zeros(size(X));
I3 = zeros(size(X));

ind1 = X > 2*b;
ind2 = and(X<=2*b,X>=b/2);
ind3 = X<b/2;
ind23 = X/2<=b;

x1 = X(ind1);
x2 = X(ind2);
x3 = X(ind3);
x23 = X(ind23);

% Case 1:

if ~isempty(x1)
    Q = fix(N*max(b./x1))+1;
    [x,w] = gaussQuadPower(alpha,Q);
    x = (b./x1).*(x');
    w = (b./x1).^(1-alpha).*w';
    fun = @(u)(1./(u-1));
    I1(ind1) = x1.^(-alpha).*sum(fun(x).*w,2);
    
end
if ~isempty(x23)
    % Common part of cases 2-3 [0,x/2]
    Q = fix(N*abs(alpha*log(min(x3))))+1;
    Q = max(Q,5);
    [x,w] = gaussQuadPower(alpha,Q);
    x = x.*1/2;
    w = w.*(1/2).^(1-alpha);
    fun1 = @(u)(1./(u-1));
    I1(ind23) = x23.^(-alpha).*sum(fun1(x).*w);
end
if ~isempty(x2)
    % case 2
    [x,w] = Gauss_Legendre1D(N,0,1);
    x = 1/2 + (b./x2 - 1/2).*x';
    w = (b./x2 - 1/2).*w';
    I2(ind2) = x2.^(-alpha).*(sum(aux(x).*w,2)+ log(2*abs(b-x2)./x2));
end
% case 3
if ~isempty(x3)
    [x,w] = Gauss_Legendre1D(Q,1/2,2);
    H = sum(aux(x).*w);
    I2(ind3) = x3.^(-alpha).*(H + log(2));
    
    % Use a log scale for the remainder
    Q = max(fix(N*(log(b/(2*min(x3)))))+1,N);
    [x,w] = Gauss_Legendre1D(Q,0,1);
    x = log(2*x3) + log(b./(2*x3)).*x';
    w = log(b./(2*x3)).*w';
    
    fun3 = @(t,u)(exp(-(alpha)*t)./(1 - exp(u - t)));
    U = log(x3).*ones(1,size(x,2));
    I3(ind3) = sum(fun3(x,U).*w,2);
end
I = I1 + I2 + I3;
%



    function[r] = aux(u)
        
        r = (u.^(-alpha) - 1)./(u-1);
        r(u==1) = -alpha;
        
    end

end



