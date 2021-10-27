function out = kernel_t2_2d(x,y,nu)
dxy = vecnorm(x,y,2);
not_singular = dxy > 1e-13;
out = not_singular .* sum((x-y).*(nu(x)-nu(y)),2)./(dxy.^2) + (~not_singular) 
end