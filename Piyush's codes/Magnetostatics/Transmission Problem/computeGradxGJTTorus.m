function out = computeGradxGJTTorus(J,R,r,Xeval)

    % Gauss legendre quadrature points
    [X,W] = lgwt(40,0,2*pi);

    % Tensorizing the quadrature for double integral 
    N = size(X,1);

    % Tensorized quadrature points and weights
    theta = repmat(X,N,1);
    phi = repelem(X,N,1);
    WW = repmat(W,N,1).*repelem(W,N,1);
    
    % Tangential source current at quadrature points
    JJ = [-sin(theta) cos(theta) 0*theta];

    % Points on the torus
    Y = [ (R+r*cos(phi)).*cos(theta) (R+r*cos(phi)).*sin(theta) r*sin(phi) ];

    % Jacobian of transformation
    surf_elem = r*(R+r*cos(phi));
    
    % Combining the source current with the jacobian
    Jtimes_surf_elem = JJ.*surf_elem;

    % Size of tensorized quadrature rule
    N = size(theta,1);

    % repeating the points X for integration.
    Neval = size(Xeval,1);
    XXeval = repelem(Xeval,N,1);
    
    % Repeating the rest of the integrand for integration
    YY = repmat(Y,Neval,1);
    Jtimes_surf_elem = repmat(Jtimes_surf_elem,Neval,1);
    gradxG = (YY - XXeval)./vecnorm(YY-XXeval,2,2).^3;

    % Integrand of size Neval * N X 3
    % has Neval repeating blocks of size N X 3 which represent the vector 
    % integrand for each evaluation point
    integrand1 = gradxG(:,1).*Jtimes_surf_elem;
    integrand2 = gradxG(:,2).*Jtimes_surf_elem;
    integrand3 = gradxG(:,3).*Jtimes_surf_elem;

    % Reshaping the integrand for integration
    integrand1 = reshape(integrand1,[N Neval*3]);
    integrand2 = reshape(integrand2,[N Neval*3]);
    integrand3 = reshape(integrand3,[N Neval*3]);

    % Evaluating all the integrals
    integral1 = sum(WW.*integrand1,1);
    integral2 = sum(WW.*integrand2,1);
    integral3 = sum(WW.*integrand3,1);

    out{1} = 1/4/pi * J * reshape(integral1,[Neval 3]);
    out{2} = 1/4/pi * J * reshape(integral2,[Neval 3]);
    out{3} = 1/4/pi * J * reshape(integral3,[Neval 3]);

end