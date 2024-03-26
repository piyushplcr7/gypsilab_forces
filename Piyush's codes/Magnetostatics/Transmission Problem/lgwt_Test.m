% script to test lgwt

[x,w] = lgwt(20,0,2*pi);

% Integrating cos^2 theta
I1 = sum(w.*cos(x).^2)

I2 = sum(w.*sin(x).^2)

[x1,w1] = lgwt(20,0,1);

I3 = sum(w1.*exp(x1))

exp(1)-1

% Gauss legendre quadrature points
    [X,W] = lgwt(20,0,2*pi);

    % Tensorizing the quadrature for double integral 
    N = size(X,1);

    theta = repmat(X,N,1);
    phi = repelem(X,N,1);

    WW = repmat(W,N,1).*repelem(W,N,1);
    
    R = 3.221;
    r = 0.799;

    area_torus = sum(WW.*(R+r*cos(phi))*r)
    analytic = 4*pi^2 * r*R





