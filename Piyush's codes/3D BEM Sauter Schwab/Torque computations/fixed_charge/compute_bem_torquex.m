function torquex = compute_bem_torquex(mesh,R,Xcg,Psi)

    S0_Gamma = fem(mesh,'P0');
    Nux = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[1 0 0],X-Xcg);

    % Kernels
    kernelx = @(x,y,z) sum(z.*(Nux(x) - Nux(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    t2matx = panel_oriented_assembly(mesh,kernelx,S0_Gamma,S0_Gamma);
    torquex = 0.5 * dot(Psi,t2matx*Psi);

end