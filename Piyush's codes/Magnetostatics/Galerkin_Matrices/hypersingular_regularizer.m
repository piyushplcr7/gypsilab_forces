function reg = hypersingular_regularizer(mesh,order)
    P0 = fem(mesh,'P0');
    P1 = fem(mesh,'P1');
    Gamma = dom(mesh,order);
    
    % Building the single layer matrix
    Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); 
    V = 1/4/pi*integral(Gamma,Gamma,P0,Gxy,P0);
    V = V + 1/4/pi*regularize(Gamma,Gamma,P0,'[1/r]',P0);

    % Solving for equilibrium density
    %M00 = integral(Gamma,P0,P0);
    M10 = integral(Gamma,P1,P0);
    N = size(M10,1);

    eq_density = V\(M10'*ones(N,1));
    reg = M10 * eq_density * eq_density' * M10';

end