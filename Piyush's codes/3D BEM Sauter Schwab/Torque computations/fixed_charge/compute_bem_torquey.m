function torquey = compute_bem_torquey(mesh,R,Xcg,Psi)

    S0_Gamma = fem(mesh,'P0');
    Nuy = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[0 1 0],X-Xcg);

    % Kernels
    kernely = @(x,y,z) sum(z.*(Nuy(x) - Nuy(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    t2maty = panel_oriented_assembly(mesh,kernely,S0_Gamma,S0_Gamma);
    torquey = 0.5 * dot(Psi,t2maty*Psi);

end