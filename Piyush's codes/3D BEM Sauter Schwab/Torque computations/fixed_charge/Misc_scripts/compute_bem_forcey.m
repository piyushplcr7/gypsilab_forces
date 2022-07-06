function forcey = compute_bem_forcey(mesh,R,Psi)

    S0_Gamma = fem(mesh,'P0');
    Nuy = @(X) (vecnorm(X,2,2)<R).* ones(size(X,1),1)*[0 1 0];

    % Kernels
    kernely = @(x,y,z) sum(z.*(Nuy(x) - Nuy(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    t2maty = panel_oriented_assembly(mesh,kernely,S0_Gamma,S0_Gamma);
    forcey = 0.5 * dot(Psi,t2maty*Psi);

end