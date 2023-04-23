% Logarithmic derivative of Gamma Function
function psi = log_derv_gamma(x,n)
    % Euler Mascheroni constant
    y = 0.57721566490153286060651209;
    
    psi = -(y+1./x);

    for i = 1:n
        psi = psi + x./(i^2+i*x);
    end

end