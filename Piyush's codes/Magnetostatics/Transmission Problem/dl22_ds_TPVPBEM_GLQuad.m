function out = dl22_ds_TPVPBEM_GLQuad(J,R,r,Xeval,Vel)

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
    xmy = XXeval - YY;
    
    gradxxG_velx = 3 * xmy .*dot(xmy,Vel(XXeval),2)./vecnorm(xmy,2,2).^5 - Vel(XXeval)./vecnorm(xmy,2,2).^3;

    % Integrand of size Neval * N X 3
    % has Neval repeating blocks of size N X 3 which represent the vector 
    % integrand for each evaluation point
    integrand = cross(gradxxG_velx,Jtimes_surf_elem,2);

    % Reshaping the integrand for integration
    integrand = reshape(integrand,[N Neval*3]);

    % Evaluating all the integrals
    integral = sum(WW.*integrand,1);
    
    % Getting columns of DHJ
    out = 1/4/pi * J * reshape(integral,[Neval 3]);

end