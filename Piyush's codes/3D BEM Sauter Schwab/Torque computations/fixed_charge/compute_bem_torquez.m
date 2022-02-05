function torquez = compute_bem_torquez(mesh,R,Xcg,Psi)

    S0_Gamma = fem(mesh,'P0');
    Nuz = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[0 0 1],X-Xcg);

    % Kernels
    kernelz = @(x,y,z) sum(z.*(Nuz(x) - Nuz(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    t2matz = panel_oriented_assembly(mesh,kernelz,S0_Gamma,S0_Gamma);
    torquez = 0.5 * dot(Psi,t2matz*Psi);

end