function torques = compute_bem_torques(mesh,R,Xcg,Psi)

    S0_Gamma = fem(mesh,'P0');

    Nux = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[1 0 0],X-Xcg);
    Nuy = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[0 1 0],X-Xcg);
    Nuz = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[0 0 1],X-Xcg);

    % Kernels
    kernelx = @(x,y,z) sum(z.*(Nux(x) - Nux(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    kernely = @(x,y,z) sum(z.*(Nuy(x) - Nuy(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    kernelz = @(x,y,z) sum(z.*(Nuz(x) - Nuz(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);

    t2matx = panel_oriented_assembly(mesh,kernelx,S0_Gamma,S0_Gamma);
    torquex = 0.5 * dot(Psi,t2matx*Psi);
    t2maty = panel_oriented_assembly(mesh,kernely,S0_Gamma,S0_Gamma);
    torquey = 0.5 * dot(Psi,t2maty*Psi);
    t2matz = panel_oriented_assembly(mesh,kernelz,S0_Gamma,S0_Gamma);
    torquez = 0.5 * dot(Psi,t2matz*Psi);

    torques = [torquex torquey torquez];

end