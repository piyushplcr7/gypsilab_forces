% Sauter and Schwab quadrature technique. Assumes that panels are given in
% the correct form. Implementation done for Triangle reference elements
function integral = sstri_integrate(kernel,hbasisx,hbasisy,chi_tau,chi_t,g_tau,g_t,relation)
    integral = 0;
    
    % definition for k3
    function val = k3(X)
        xh = [X(1);X(2)];
        yh = [X(3);X(4)];
        val = hbasisx(xh) * hbasisy(yh) * kernel(chi_tau(xh),chi_t(yh),chi_t(yh)-chi_tau(xh)) * g_tau(xh) * g_t(yh);
    end

    % Definition of the integrand for identical panels case
    function val = integrand1(n1,n2,n3,z)
        
        val = k3( z * [1;
            1-n1+n1*n2;
            1-n1*n2*n3;
            1-n1] );
        
        val = val + k3( z * [1-n1*n2*n3;
            1-n1;
            1;
            1-n1+n1*n2] );
        
        val = val + k3( z * [1;
            n1*(1-n2+n2*n3);
            1-n1*n2;
            n1*(1-n2)] );
        
        val = val + k3( z * [1-n1*n2;
            n1*(1-n2);
            1;
            n1*(1-n2+n2*n3)] );
        
        val = val + k3( z * [1-n1*n2*n3;
            n1*(1-n2*n3);
            1;
            n1*(1-n2)] );
        
        val = val + k3( z * [1;
            n1*(1-n2);
            1-n1*n2*n3;
            n1*(1-n2*n3)] );
        
        val = val * z^3 * n1^2 * n2;
    end
    
    % Integrand for common edge case
    function val = integrand2(n1,n2,n3,z)
        
        val = k3( z * [1;
            n1;
            1-n1*n2*n3;
            n1*n2*(1-n3)] );
        
        val = val + k3( z * [1-n1*n2;
            n1*(1-n2);
            1;
            n1*n2*n3] );
        
        val = val + k3( z * [1-n1*n2*n3;
            n1*n2*(1-n3);
            1;
            n1] );
        
        val = val + k3( z * [1-n1*n2*n3;
            n1*(1-n2*n3);
            1;
            n1*n2] );
        
        val = val * z^3 * n1^2 * n2;
        
        val = val + z^3 * n1^2 * k3( z * [1;
            n1*n3;
            1-n1*n2;
            n1*(1-n2)] );
    end

    % Integrand for common vertex case
    function val = integrand3(n1,n2,n3,z)
        
        val = k3( z * [1;
            n1;
            n2;
            n2*n3] );
        
        val = val + k3( z * [n2;
            n2*n3;
            1;
            n1] );
        
        val = val * z^3 * n2;
    end

% Integrand for far away panel case
    function val = integrand4(n1,n2,n3,n4)
        
        val = k3([n1;
            n2 * n1;
            n3;
            n4 * n3] );
        
        val = val * n1 * n3;
    end
    
    switch relation
        case "identical"
            disp("identical panels case");
            integral = integrate4dim(@integrand1);
        case "common_edge"
            disp("common edge case");
            integral = integrate4dim(@integrand2);
        case "common_vertex"
            disp("common vertex case");
            integral = integrate4dim(@integrand3);
        case "far_away"
            disp("far away panels case");     
            integral = integrate4dim(@integrand4);
        otherwise
            disp("Unrecognized case, terminating without computation");
    end
    
end