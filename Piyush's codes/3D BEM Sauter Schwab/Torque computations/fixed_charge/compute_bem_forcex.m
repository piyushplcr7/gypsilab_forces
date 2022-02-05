function forcex = compute_bem_forcex(mesh,R,Psi)

    S0_Gamma = fem(mesh,'P0');

    Nux = @(X) (vecnorm(X,2,2)<R).* ones(size(X,1),1)*[1 0 0];

    % Kernels
    kernelx = @(x,y,z) sum(z.*(Nux(x) - Nux(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    t2matx = panel_oriented_assembly(mesh,kernelx,S0_Gamma,S0_Gamma);
    forcex = 0.5 * dot(Psi,t2matx*Psi);

end