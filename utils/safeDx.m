function [r] = safeDx(f,I,dx,x)
% Compute approximate derivative while staying in interval of definition

if nargin == 2
    r = @(x)(safeDx(f,I,[],x));
elseif nargin == 3
    r = @(x)(safeDx(f,I,dx,x));
else
    a = I(1); b = I(2);
    
    if ~exist('dx','var')||isempty(dx)
        dx = min(1e-5*(b-a),1e-5);
    end
    
    
    R = @(x)((f(x + (b-x)/(b-a)*dx) - f(x - (x - a)/(b-a)*dx))./dx);
    t = linspace(a,b,1000);
    F = fft(R(t));
    r = (f(x + (b-x)/(b-a)*dx) - f(x - (x - a)/(b-a)*dx))./dx;

    
    
end




end

