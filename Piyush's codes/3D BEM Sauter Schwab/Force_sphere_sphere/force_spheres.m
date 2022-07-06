% Force formula for two spheres at a constant potential difference enforced
% by a battery
function out = force_spheres(d,R,epsilon,V0)
    alpha0 = acosh(d/(2*R));
    out = pi/4 * epsilon * V0^2 * Gsum(300,alpha0);
end