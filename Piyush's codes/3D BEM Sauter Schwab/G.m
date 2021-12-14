function out = G(n,alpha)
    out = exp(-(n+0.5).*alpha)./sinh((n+0.5)*alpha);
end