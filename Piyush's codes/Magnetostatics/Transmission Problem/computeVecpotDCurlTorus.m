function DHj = computeVecpotDCurlTorus(J,R,r,Xeval)

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
    del1gradxG = 3* xmy .* xmy(:,1)./vecnorm(xmy,2,2).^5 - [1 0 0]./vecnorm(xmy,2,2).^3;
    del2gradxG = 3* xmy .* xmy(:,2)./vecnorm(xmy,2,2).^5 - [0 1 0]./vecnorm(xmy,2,2).^3;
    del3gradxG = 3* xmy .* xmy(:,3)./vecnorm(xmy,2,2).^5 - [0 0 1]./vecnorm(xmy,2,2).^3;

    % Integrand of size Neval * N X 3
    % has Neval repeating blocks of size N X 3 which represent the vector 
    % integrand for each evaluation point
    integrand1 = cross(del1gradxG,Jtimes_surf_elem,2);
    integrand2 = cross(del2gradxG,Jtimes_surf_elem,2);
    integrand3 = cross(del3gradxG,Jtimes_surf_elem,2);

    % Reshaping the integrand for integration
    integrand1 = reshape(integrand1,[N Neval*3]);
    integrand2 = reshape(integrand2,[N Neval*3]);
    integrand3 = reshape(integrand3,[N Neval*3]);

    % Evaluating all the integrals
    integral1 = sum(WW.*integrand1,1);
    integral2 = sum(WW.*integrand2,1);
    integral3 = sum(WW.*integrand3,1);
    
    % Getting columns of DHJ
    DHj1 = 1/4/pi * J * reshape(integral1,[Neval 3]);
    DHj2 = 1/4/pi * J * reshape(integral2,[Neval 3]);
    DHj3 = 1/4/pi * J * reshape(integral3,[Neval 3]);

    DHj{1} = [DHj1(:,1) DHj2(:,1) DHj3(:,1)];
    DHj{2} = [DHj1(:,2) DHj2(:,2) DHj3(:,2)];
    DHj{3} = [DHj1(:,3) DHj2(:,3) DHj3(:,3)];

end