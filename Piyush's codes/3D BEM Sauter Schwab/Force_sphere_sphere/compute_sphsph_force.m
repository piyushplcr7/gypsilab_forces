% Force formula for two identical spheres with one connected to ground
function force = compute_sphsph_force(r,s,Q)
    % Euler Mascheroni constant
    y = .57721566490153286060651209;
    psi_at_half =  -2*log(2) -y;
    force = -1/(2*r) * Q^2 * 1/s * 1/(0.5*log(r/s) - psi_at_half)^2;

end