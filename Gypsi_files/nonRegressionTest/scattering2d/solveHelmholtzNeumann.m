function [mu] = solveHelmholtzNeumann(Gamma,Vh,dudn,k)


G = @(X,Y)(femGreenKernel(X,Y,'[H0(kr)]',k));

H = 1i/4 .* (k^2 * integral(Gamma,Gamma,ntimes(Vh),G,ntimes(Vh)) ...
    - integral(Gamma,Gamma,nxgrad(Vh),G,nxgrad(Vh)));

Hr = -1/(2*pi) .*  (k^2 * regularize(Gamma,Gamma,ntimes(Vh),'[log(r)]',ntimes(Vh)) ...
    - regularize(Gamma,Gamma,nxgrad(Vh),'[log(r)]',nxgrad(Vh)));

H = H+Hr;

mu = H \dudn;

end

