function torques = compute_bem_forces_gypsi(mesh,R,Psi)

    S0_Gamma = fem(mesh,'P0');

    Nux = @(X) (vecnorm(X,2,2)<R).* ones(size(X,1),1)*[1 0 0];
    Nuy = @(X) (vecnorm(X,2,2)<R).* ones(size(X,1),1)*[0 1 0];
    Nuz = @(X) (vecnorm(X,2,2)<R).* ones(size(X,1),1)*[0 0 1];

    % Kernels
    kernelx = @(x,y,z) sum(z.*(Nux(x) - Nux(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    kernely = @(x,y,z) sum(z.*(Nuy(x) - Nuy(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    kernelz = @(x,y,z) sum(z.*(Nuz(x) - Nuz(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);

    t2matx = panel_oriented_assembly(mesh,kernelx,S0_Gamma,S0_Gamma);
    forcex = 0.5 * dot(Psi,t2matx*Psi);
    t2maty = panel_oriented_assembly(mesh,kernely,S0_Gamma,S0_Gamma);
    forcey = 0.5 * dot(Psi,t2maty*Psi);
    t2matz = panel_oriented_assembly(mesh,kernelz,S0_Gamma,S0_Gamma);
    forcez = 0.5 * dot(Psi,t2matz*Psi);

    torques = [forcex forcey forcez];

end