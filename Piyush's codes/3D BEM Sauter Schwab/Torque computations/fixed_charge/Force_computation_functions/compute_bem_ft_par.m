function [fbem,tbem] = compute_bem_ft_par(mesh,R,Xcg,Psi)

    S0_Gamma = fem(mesh,'P0');

    % Translation fields
    Nux = @(X) (vecnorm(X,2,2)<R).* ones(size(X,1),1)*[1 0 0];
    Nuy = @(X) (vecnorm(X,2,2)<R).* ones(size(X,1),1)*[0 1 0];
    Nuz = @(X) (vecnorm(X,2,2)<R).* ones(size(X,1),1)*[0 0 1];

    % Rotational fields
    Nuxr = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[1 0 0],X-Xcg);
    Nuyr = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[0 1 0],X-Xcg);
    Nuzr = @(X) (vecnorm(X,2,2)<R).* cross(ones(size(X,1),1)*[0 0 1],X-Xcg);

    % Kernels translation
    kernelx = @(x,y,z) sum(z.*(Nux(x) - Nux(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    kernely = @(x,y,z) sum(z.*(Nuy(x) - Nuy(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    kernelz = @(x,y,z) sum(z.*(Nuz(x) - Nuz(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);

    % Kernels rotation
    kernelxr = @(x,y,z) sum(z.*(Nuxr(x) - Nuxr(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    kernelyr = @(x,y,z) sum(z.*(Nuyr(x) - Nuyr(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);
    kernelzr = @(x,y,z) sum(z.*(Nuzr(x) - Nuzr(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);

    kernels = cell(6,1);
    kernels{1} = kernelx; 
    kernels{2} = kernely;
    kernels{3} = kernelz;
    kernels{4} = kernelxr;
    kernels{5} = kernelyr;
    kernels{6} = kernelzr;
    ft = cell(6,1);

    parfor  i = 1:6
        t2mat = panel_oriented_assembly(mesh,kernels{i},S0_Gamma,S0_Gamma);
        ft{i} = 0.5 * dot(Psi,t2mat*Psi);
    end

    fbem = [ft{1} ft{2} ft{3}];
    tbem = [ft{4} ft{5} ft{6}];

end