% Function to evaluate the limit encountered while considering field energy

function out = limitfunction(R,gamma)
    
alpha = R.^2+gamma.^2;
beta = 2 * gamma.*R./alpha;

out = -4/3 * pi * R.^2./sqrt(alpha) .* (sqrt(1+beta).*(beta-2) + sqrt(1-beta).*(beta+2))./beta.^2;

end